"""
t3 utils writer module
"""

import os
import shutil

from mako.template import Template

from t3.common import get_rmg_species_from_a_species_dict
from t3.utils.generator import generate_radicals

METHOD_MAP = {'CSE': 'chemically-significant eigenvalues',
              'RS': 'reservoir state',
              'MSC': 'modified strong collision',
              }


def write_rmg_input_file(rmg: dict,
                         t3: dict,
                         iteration: int,
                         path: str,
                         walltime: str = '00:00:00:00',
                         ):
    """
    Write an RMG input file to the given file path.
    Will create the directory if needed.

    Args:
        rmg (dict): The arguments to write in a keyword argument dictionary format.
        t3 (dict): The T3 arguments in a keyword argument dictionary format. Includes atol and rtol for SA.
        iteration (int): The T3 iteration, used to determine ``core_tolerance`` and ``tolerance_interrupt_simulation``.
        path (str): The path where the RMG input file should be saved.
        walltime (str, optional): The time cap for an RMG run. Should pass here t3['options']['max_RMG_walltime']
    """
    rmg = rmg.copy()
    rmg_input = ''
    iteration -= 1  # iteration is 1-indexed, convert to 0-indexed for list indexing

    # database
    database = rmg['database']
    # the following args type could be either str or list, detect str and format accordingly
    if isinstance(database['kinetics_depositories'], str) and database['kinetics_depositories'][0] != "'":
        database['kinetics_depositories'] = f"'{database['kinetics_depositories']}'"
    if isinstance(database['kinetics_estimator'], str) and database['kinetics_estimator'][0] != "'":
        database['kinetics_estimator'] = f"'{database['kinetics_estimator']}'"
    if isinstance(database['kinetics_families'], str) and database['kinetics_families'][0] != "'":
        database['kinetics_families'] = f"'{database['kinetics_families']}'"
    database_template = """database(
    thermoLibraries=${thermo_libraries},
    reactionLibraries=${kinetics_libraries},
    transportLibraries=${transport_libraries},
    seedMechanisms=${seed_mechanisms},
    kineticsDepositories=${kinetics_depositories},
    kineticsFamilies=${kinetics_families},
    kineticsEstimator=${kinetics_estimator},
)
"""
    rmg_input += Template(database_template).render(**database)

    # species
    species = rmg['species']
    species_template = """
species(
    label='${label}',
    reactive=${reactive},
    structure=${structure},
)
"""
    for spc in species:
        if spc['adjlist'] is not None:
            structure = f'adjacencyList("""'
            structure += spc['adjlist']
            structure += '""")'
        elif spc['smiles'] is not None:
            structure = f"SMILES('{spc['smiles']}')"
        elif spc['inchi'] is not None:
            structure = f"InChI('{spc['inchi']}')"
        else:
            raise ValueError(f"A species must have either an adjlist, smiles, or inchi descriptor.")
        rmg_input += Template(species_template).render(label=spc['label'],
                                                       reactive=spc['reactive'],
                                                       structure=structure)
        if spc['seed_all_rads'] is not None:
            species_to_process = get_rmg_species_from_a_species_dict(species_dict=spc, raise_error=False)
            if species_to_process is not None:
                radical_tuples = generate_radicals(species=species_to_process,
                                                   types=spc['seed_all_rads'],
                                                   )
                for radical_tuple in radical_tuples:
                    rmg_input += Template(species_template).render(label=radical_tuple[0],
                                                                   reactive=True,
                                                                   structure=f"SMILES('{radical_tuple[1]}')")

    # reactors
    reactors = rmg['reactors']
    gas_batch_constant_t_p_template_template = """
simpleReactor(
    temperature=${temperature},
    pressure=${pressure},
    initialMoleFractions={${concentrations()}    },
    ${termination}
    nSims=${conditions_per_iteration},${balance}${constant}
)
<%def name="concentrations()">
% for spc in species_list:
        '${spc["label"]}': ${spc["concentration"]},
% endfor
</%def>
"""
    liquid_batch_constant_t_v_template = """
liquidReactor(
    temperature=${temperature},
    initialConcentrations={${concentrations()}    },
    ${termination}
    nSims=${conditions_per_iteration},${constant}
)
<%def name="concentrations()">
% for spc in species_list:
        '${spc["label"]}': (${spc["concentration"]}, 'mol/cm^3'),
% endfor
</%def>
"""
    for reactor in reactors:
        if isinstance(reactor['T'], float):
            temperature = f"({reactor['T']}, 'K')"
        elif isinstance(reactor['T'], list):
            temperature = [(t, 'K') for t in reactor['T']]
        else:
            raise ValueError(f"The reactor temperature must be a float or a list,\n"
                             f"got {reactor['T']} which is a {type(reactor['T'])}.")
        species_list = [{'label': spc['label'], 'concentration': spc['concentration']} for spc in species]
        species_list.sort(key=lambda spc: spc['concentration'], reverse=True)
        termination = ''
        if reactor['termination_conversion'] is not None:
            termination += f"terminationConversion={reactor['termination_conversion']},"
        if reactor['termination_time'] is not None:
            termination += '\n    ' if termination else ''
            termination += f"terminationTime={reactor['termination_time']},"
        if reactor['termination_rate_ratio'] is not None:
            termination += '\n    ' if termination else ''
            termination += f"terminationRateRatio={reactor['termination_rate_ratio']},"
        constant = ''
        for spc in species:
            if spc['constant']:
                if not constant:
                    constant = '\n    constantSpecies=['
                constant += f"'{spc['label']}', "
        constant += '],' if constant else ''

        if reactor['type'] == 'gas batch constant T P':
            if isinstance(reactor['P'], float):
                pressure = f"({reactor['P']}, 'bar')"
            elif isinstance(reactor['P'], list):
                pressure = [(p, 'K') for p in reactor['P']]
            else:
                raise ValueError(f"The reactor pressure must be a float or a list,\n"
                                 f"got {reactor['P']} which is a {type(reactor['P'])}.")
            balance = ''
            for spc in species:
                if spc['balance']:
                    balance = f"\n    balanceSpecies='{spc['label']}',"
                    break
            rmg_input += Template(gas_batch_constant_t_p_template_template).render(
                temperature=temperature,
                pressure=pressure,
                species_list=species_list,
                termination=termination,
                conditions_per_iteration=reactor['conditions_per_iteration'],
                balance=balance,
                constant=constant,
            )

        elif reactor['type'] == 'liquid batch constant T V':
            rmg_input += Template(liquid_batch_constant_t_v_template).render(
                temperature=temperature,
                species_list=species_list,
                termination=termination,
                conditions_per_iteration=reactor['conditions_per_iteration'],
                constant=constant,
            )

    # solvent
    solvent_template = """solvation(solvent='${solvent}')

"""
    solvent = ''
    for spc in species:
        # the schema assures that there's only one species defined as the solvent
        # TODO: assure that the requested solvent actually exists in the RMG database
        if spc['solvent']:
            solvent = spc['label']
            break

    if solvent:
        rmg_input += Template(solvent_template).render(
            solvent=solvent,
        )

    # model
    model_input = rmg['model']
    model_template = """model(
    toleranceMoveToCore=${tol_move_to_core},
    toleranceInterruptSimulation=${tolerance_interrupt_simulation},${args}
)
"""
    model = dict()
    model['tol_move_to_core'] = model_input['core_tolerance'][iteration] \
        if len(model_input['core_tolerance']) >= iteration + 1 else model_input['core_tolerance'][-1]
    model['tolerance_interrupt_simulation'] = model_input['tolerance_interrupt_simulation'][iteration] \
        if len(model_input['tolerance_interrupt_simulation']) >= iteration + 1 \
        else model_input['tolerance_interrupt_simulation'][-1]
    model_keys_to_skip = ['core_tolerance', 'tolerance_interrupt_simulation', 'atol', 'rtol', 'sens_atol', 'sens_rtol']
    args = ''
    for key, value in model_input.items():
        if key not in model_keys_to_skip and value is not None:
            args += f"\n    {to_camel_case(uv=key)}={value},"
    model['args'] = args
    rmg_input += Template(model_template).render(**model)

    # simulator
    if t3['sensitivity'] is not None and t3['sensitivity']['adapter'] == 'RMG':
        simulator_template = """\nsimulator(atol=${atol}, rtol=${rtol}, sens_atol=${sens_atol}, sens_rtol=${sens_rtol})\n"""
        rmg_input += Template(simulator_template).render(atol=model_input['atol'],
                                                         rtol=model_input['rtol'],
                                                         sens_atol=t3['sensitivity']['atol'],
                                                         sens_rtol=t3['sensitivity']['rtol']
                                                         )
    else:
        simulator_template = """\nsimulator(atol=${atol}, rtol=${rtol})\n"""
        rmg_input += Template(simulator_template).render(atol=model_input['atol'],
                                                         rtol=model_input['rtol'],
                                                         )

    # pressureDependence
    if rmg['pdep'] is not None:
        pdep = rmg['pdep'].copy()
        pdep_template = """
pressureDependence(
    method='${method}',
    maximumGrainSize=(${max_grain_size}, 'kJ/mol'),
    minimumNumberOfGrains=${max_number_of_grains},
    temperatures=(${T_min}, ${T_max}, 'K', ${T_count}),
    pressures=(${P_min}, ${P_max}, 'bar', ${P_count}),
    interpolation=${interpolation},
    maximumAtoms=${max_atoms},
)
"""
        pdep['method'] = METHOD_MAP[pdep['method']] if pdep['method'] not in METHOD_MAP.values() else pdep['method']
        pdep['T_min'], pdep['T_max'], pdep['T_count'] = pdep['T']
        pdep['P_min'], pdep['P_max'], pdep['P_count'] = pdep['P']
        del pdep['T']
        del pdep['P']
        if pdep['interpolation'] == 'PDepArrhenius':
            pdep['interpolation'] = ('PDepArrhenius',)
        else:
            pdep['interpolation'] = ('Chebyshev', pdep['T_basis_set'], pdep['P_basis_set'])
        del pdep['T_basis_set']
        del pdep['P_basis_set']
        rmg_input += Template(pdep_template).render(**pdep)

    # options
    options = rmg['options']
    if options is not None:
        options_template = """
options(
    name='${seed_name}',
    generateSeedEachIteration=${generate_seed_each_iteration},
    saveSeedToDatabase=${save_seed_to_database},
    units='${units}',
    generateOutputHTML=${save_html},
    generatePlots=${generate_plots},
    saveSimulationProfiles=${save_simulation_profiles},
    verboseComments=${verbose_comments},
    saveEdgeSpecies=${save_edge},
    keepIrreversible=${keep_irreversible},
    trimolecularProductReversible=${trimolecular_product_reversible},
    wallTime='${walltime}',
    saveSeedModulus=${save_seed_modulus},
)
"""
        options['walltime'] = walltime
        rmg_input += Template(options_template).render(**options)

    # generatedSpeciesConstraints
    species_constraints = rmg['species_constraints']
    if species_constraints is not None:
        species_constraints_template = """
generatedSpeciesConstraints(
    allowed=${allowed},
    maximumCarbonAtoms=${max_C_atoms},
    maximumOxygenAtoms=${max_O_atoms},
    maximumNitrogenAtoms=${max_N_atoms},
    maximumSiliconAtoms=${max_Si_atoms},
    maximumSulfurAtoms=${max_S_atoms},
    maximumHeavyAtoms=${max_heavy_atoms},
    maximumRadicalElectrons=${max_radical_electrons},
    maximumSingletCarbenes=${max_singlet_carbenes},
    maximumCarbeneRadicals=${max_carbene_radicals},
    allowSingletO2=${allow_singlet_O2},
)
"""
        rmg_input += Template(species_constraints_template).render(**species_constraints)

    if not os.path.isdir(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    with open(path, 'w') as f:
        f.writelines(rmg_input)


def write_pdep_network_file(network_name: str,
                            method: str,
                            pdep_sa_path: str,
                            rmg_pdep_path: str,
                            ) -> tuple:
    """
    Adding a P-dep SA directive to an Arkane network input file.

    Args:

        network_name (str): The name of the original network file, e.g. 'network32_1'.
        method (str): 'CSE', 'MSC' or 'RS'.
        pdep_sa_path (str): The path to the PDep_SA iteration folder.
        rmg_pdep_path (str): The path to the RMG/pdep iteration folder.

    Raises:
        ValueError: If T/P ranges could not be parsed from the file.

    Returns:
        tuple: Isomer labels of the current network.
    """
    # copy the network file into a new folder (and rename it to input.py)
    sa_pdep_path = os.path.join(pdep_sa_path, network_name, method)
    if not os.path.isdir(sa_pdep_path):
        os.makedirs(sa_pdep_path)
    input_file_path = os.path.join(sa_pdep_path, 'input.py')
    shutil.copyfile(src=os.path.join(rmg_pdep_path, network_name + '.py'),
                    dst=input_file_path)

    with open(input_file_path, 'r') as f:
        lines = f.readlines()
        new_lines, isomer_labels = list(), list()
        t_min, t_max, p_min, p_max = None, None, None, None
        parse_tp, parse_isomers = False, (False, False)
        for line in lines:
            if 'pressureDependence(' in line:
                parse_tp = True
            if 'network(' in line:
                parse_isomers = (True, False)
            if parse_isomers[0] and 'isomers =' in line:
                parse_isomers = (True, True)
            if 'reactants =' in line:
                parse_isomers = (False, False)
            if parse_tp:
                if 'Tmin' in line:
                    #     Tmin = (300,'K'),
                    t_min = line.split('(')[1].split(',')[0]
                elif 'Tmax' in line:
                    #     Tmax = (2200,'K'),
                    t_max = line.split('(')[1].split(',')[0]
                elif 'Pmin' in line:
                    #     Pmin = (0.01,'bar'),
                    p_min = line.split('(')[1].split(',')[0]
                elif 'Pmax' in line:
                    #     Pmax = (100,'bar'),
                    p_max = line.split('(')[1].split(',')[0]
            if all(parse_isomers) and "'," in line:
                #         'C=O(26)',
                isomer_labels.append(line.split("'")[1])
            if 'method = ' in line:
                #     method = 'chemically-significant eigenvalues',
                splits = line.split("'")
                new_lines.append(f"{splits[0]}'{METHOD_MAP[method]}'{splits[2]}")
            elif 'rmgmode' in line:
                new_lines.append(line)
                if any(param is None for param in [t_min, t_max, p_min, p_max]):
                    raise ValueError(f'Could not parse all T/P parameters, got:\n'
                                     f'T min = {t_min}, T max = {t_max}, P min = {p_min}, P max = {p_max}.')
                sa_conditions = f"""    sensitivity_conditions = [[({t_min}, 'K'), ({p_min}, 'bar')],
                              [({t_max}, 'K'), ({p_min}, 'bar')],
                              [({t_min}, 'K'), ({p_max}, 'bar')],
                              [({t_max}, 'K'), ({p_max}, 'bar')]],"""
                new_lines.append(sa_conditions)
            else:
                new_lines.append(line)

    with open(input_file_path, 'w') as f:
        f.writelines(new_lines)

    return tuple(isomer_labels)


def to_camel_case(uv: str) -> str:
    """
    Convert an underscore variable to a camel case variable

    Args:
        uv: The underscore variable

    Returns:
        str: The camel case variable.
    """
    ccv = ''
    capitalize = False
    for char in uv:
        if char != '_':
            ccv += char.capitalize() if capitalize else char
            capitalize = False
        else:
            capitalize = True
    return ccv
