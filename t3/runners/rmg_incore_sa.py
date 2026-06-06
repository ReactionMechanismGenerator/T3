#!/usr/bin/env python3
# encoding: utf-8

"""
A script for running RMG Sensitivity Analysis (SA) incore using the rmg_env.
"""

import argparse
import os
import sys
import traceback
import yaml
import pandas as pd
import shutil
from typing import Dict, Any

try:
    from rmgpy.rmg.main import RMG
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.chemkin import load_chemkin_file
    from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
except ImportError:
    sys.stderr.write("Error: Could not import rmgpy. Ensure this script runs in the rmg_env.\n")
    sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Run RMG Sensitivity Analysis')
    parser.add_argument('-i', '--input', type=str, required=True, help='RMG input file')
    parser.add_argument('-c', '--chemkin', type=str, required=True, help='Chemkin file')
    parser.add_argument('-d', '--dictionary', type=str, required=True, help='Species dictionary')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output YAML file')
    parser.add_argument('-obs', '--observables', type=str, nargs='+', default=[], help='Observable species labels')
    parser.add_argument('-t', '--threshold', type=float, default=1e-3, help='Sensitivity threshold')
    return parser.parse_args()


def process_sa_csvs(output_dir: str) -> Dict[str, Any]:
    """Parse output CSVs to YAML dict."""
    solver_path = os.path.join(output_dir, 'solver')
    if not os.path.exists(solver_path):
        # If no solver dir, return empty (implies no SA performed)
        return {'kinetics': {}, 'thermo': {}, 'time': []}

    sa_files = [f for f in os.listdir(solver_path) if 'sensitivity' in f and f.endswith(".csv")]
    sa_dict = {'kinetics': {}, 'thermo': {}, 'time': []}
    time_set = False

    for sa_file in sa_files:
        df = pd.read_csv(os.path.join(solver_path, sa_file))
        if not time_set:
            for col in df.columns:
                if 'Time' in col:
                    sa_dict['time'] = df[col].values.tolist()
                    time_set = True
                    break

        for header in df.columns:
            if 'Time' in header: continue

            # Header format: dln[Label(ID)]/dln[k#] or dln[Label(ID)]/dG[Label(ID)]
            sa_type = None
            if '/dln[k' in header:
                sa_type = 'kinetics'
            elif '/dG[' in header:
                sa_type = 'thermo'

            if sa_type:
                # Parse Observable Label: "dln[Name(ID)]..." -> "Name(ID)"
                try:
                    # Robust split: Find string between first '[' and first ']'
                    obs_part = header.split('[', 1)[1].split(']', 1)[0]
                    # T3 expects just the label. RMG adds '(ID)'.
                    # If we use_chemkin_names=True, the label should be clean from chemkin.
                    # But RMG often appends IDs in output headers.
                    # We will strip the (ID) suffix if present to match T3 expectations.
                    if '(' in obs_part and obs_part.endswith(')'):
                        observable_label = obs_part.rsplit('(', 1)[0]
                    else:
                        observable_label = obs_part
                except IndexError:
                    observable_label = header

                if observable_label not in sa_dict[sa_type]:
                    sa_dict[sa_type][observable_label] = {}

                # Parse Parameter: "dln[k10]" -> 10, "dG[Name]" -> Name
                try:
                    # Get the part after the last '/'
                    denom = header.rsplit('/', 1)[1]
                    # Get content inside brackets
                    param_raw = denom.split('[', 1)[1].split(']', 1)[0]

                    if sa_type == 'kinetics' and param_raw.startswith('k'):
                        parameter = int(param_raw[1:])  # k10 -> 10
                    else:
                        parameter = param_raw
                except (IndexError, ValueError):
                    parameter = header

                sa_dict[sa_type][observable_label][parameter] = df[header].values.tolist()

    return sa_dict


def main():
    args = parse_arguments()

    try:
        # 1. Initialize & Load Input
        rmg = RMG(input_file=args.input, output_directory=os.path.dirname(args.input))
        rmg.load_input(args.input)

        # 2. Load Chemkin Model (Official Species)
        species_list, reaction_list = load_chemkin_file(
            args.chemkin,
            args.dictionary,
            check_duplicates=False,
            use_chemkin_names=True
        )

        # Setup RMG object structure
        rmg.reaction_model = CoreEdgeReactionModel()
        rmg.reaction_model.core.species = species_list
        rmg.reaction_model.core.reactions = reaction_list

        # 3. Apply Sensitivity Settings (CLI Override)
        # Create a map for easy lookup
        species_map = {s.label: s for s in species_list}

        sensitive_species = []
        if args.observables:
            print(f"Applying sensitivity for observables: {args.observables}")
            for obs in args.observables:
                if obs in species_map:
                    sensitive_species.append(species_map[obs])
                else:
                    print(f"Warning: Observable '{obs}' not found in Chemkin model.")

        # 4. Initialize & Update Reactors
        for reactor in rmg.reaction_systems:
            # Reconcile Initial Conditions
            new_initials = {}
            for spc, val in reactor.initial_mole_fractions.items():
                if spc.label in species_map:
                    new_initials[species_map[spc.label]] = val
                else:
                    new_initials[spc] = val
            reactor.initial_mole_fractions = new_initials

            # Force Sensitivity from CLI if provided
            if sensitive_species:
                reactor.sensitive_species = sensitive_species
                reactor.sensitivity_threshold = args.threshold
            # Else: leave whatever was in input file (but it might be broken/unreconciled objects)

            # Initialize
            settings = rmg.simulator_settings_list[-1] if rmg.simulator_settings_list else SimulatorSettings()
            reactor.initialize_model(
                core_species=rmg.reaction_model.core.species,
                core_reactions=rmg.reaction_model.core.reactions,
                edge_species=[],
                edge_reactions=[],
                pdep_networks=[],
                atol=settings.atol,
                rtol=settings.rtol
            )

        # 5. Simulate
        solver_dir = os.path.join(rmg.output_directory, 'solver')
        if os.path.exists(solver_dir): shutil.rmtree(solver_dir)
        os.makedirs(solver_dir)

        sim_settings = rmg.simulator_settings_list[-1] if rmg.simulator_settings_list else SimulatorSettings()
        mod_settings = rmg.model_settings_list[-1] if rmg.model_settings_list else ModelSettings()

        count = 0
        for index, reactor in enumerate(rmg.reaction_systems):
            if reactor.sensitive_species:
                count += 1
                # Define worksheet paths
                ws = [os.path.join(solver_dir, f'sensitivity_{index + 1}_SPC_{s.index}.csv')
                      for s in reactor.sensitive_species]

                reactor.simulate(
                    core_species=rmg.reaction_model.core.species,
                    core_reactions=rmg.reaction_model.core.reactions,
                    edge_species=[], edge_reactions=[], surface_species=[], surface_reactions=[], pdep_networks=[],
                    sensitivity=True,
                    sens_worksheet=ws,
                    model_settings=mod_settings,
                    simulator_settings=sim_settings,
                    conditions=reactor.sens_conditions if hasattr(reactor, 'sens_conditions') else None
                )

        if count == 0:
            print("No reactors had sensitive species configured.")

        # 6. Output
        sa_results = process_sa_csvs(rmg.output_directory)
        with open(args.output, 'w') as f:
            yaml.dump(sa_results, f)

    except Exception as e:
        sys.stderr.write(f"RMG SA Failed: {e}\n")
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
