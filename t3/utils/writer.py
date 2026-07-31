"""
t3 utils writer module
"""

import ast
import logging
import os
import re
import shutil
from dataclasses import dataclass

from mako.template import Template

from arc.species.perceive import perceive_molecule_from_xyz

from t3.chem import T3Species
from t3.utils.generator import generate_radicals
from t3.utils.network_thermo import TGridClampRecord, format_skipped_species, network_thermo_t_max

METHOD_MAP = {'CSE': 'chemically-significant eigenvalues',
              'RS': 'reservoir state',
              'MSC': 'modified strong collision',
              }

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ArkaneNetworkWriteResult:
    """
    The result of writing an Arkane network input file: the isomer labels the writer's own
    callers already needed, plus the T-grid clamp provenance that used to only exist as a
    ``logger.warning(...)`` line and did not survive past the run. Bundled together (rather than
    changing ``write_arkane_network_input_file`` to return a bare tuple of two things) so that
    every existing caller's ``x = write_arkane_network_input_file(...)`` call site fails loudly
    and unambiguously if it is not updated to use ``.isomer_labels``, instead of silently
    unpacking a differently-shaped tuple.

    Attributes:
        isomer_labels (tuple): Isomer labels of the current network, as previously returned bare.
        t_grid_clamp (TGridClampRecord): Provenance for whether the T grid written to
            ``dest_path`` was clamped down from what was requested. See ``TGridClampRecord``'s
            docstring for the full three-state ("clamped" / "not clamped" / caller never asked,
            i.e. this whole result is simply absent) design rationale.
    """
    isomer_labels: tuple
    t_grid_clamp: TGridClampRecord


def format_clamped_t_max(value: float) -> str:
    """
    Format a clamped Tmax value for injection into an Arkane network input file's text.

    Matches the plain-integer style RMG's own writer uses for a P-dep network's Tmax (e.g.
    ``3000`` rather than ``3000.0``) whenever the clamped value is a whole number, so a clamp does
    not introduce a stylistic difference from the numbers RMG itself would have written; a
    fractional species thermo ceiling (e.g. ``957.493``) is preserved exactly via ``repr``.

    Args:
        value (float): The clamped Tmax value, in Kelvin.

    Returns:
        str: The formatted value.
    """
    if value == int(value):
        return str(int(value))
    return repr(value)


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
    if isinstance(database['kinetics_depositories'], str) and database['kinetics_depositories'] and database['kinetics_depositories'][0] != "'":
        database['kinetics_depositories'] = f"'{database['kinetics_depositories']}'"
    if isinstance(database['kinetics_estimator'], str) and database['kinetics_estimator'] and database['kinetics_estimator'][0] != "'":
        database['kinetics_estimator'] = f"'{database['kinetics_estimator']}'"
    if isinstance(database['kinetics_families'], str) and database['kinetics_families'] and database['kinetics_families'][0] != "'":
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
            structure = 'adjacencyList("""'
            structure += spc['adjlist']
            structure += '""")'
        elif spc['smiles'] is not None:
            structure = f"SMILES('{spc['smiles']}')"
        elif spc['inchi'] is not None:
            structure = f"InChI('{spc['inchi']}')"
        else:
            raise ValueError("A species must have either an adjlist, smiles, or inchi descriptor.")
        rmg_input += Template(species_template).render(label=spc['label'],
                                                       reactive=spc['reactive'],
                                                       structure=structure)
        if spc['seed_all_rads'] is not None:
            species_to_process = get_species_obj_from_a_species_dict(species_dict=spc, raise_error=False)
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
    gas_batch_constant_t_p_template = """
simpleReactor(
    temperature=${temperature},
    pressure=${pressure},
    initialMoleFractions={${concentrations()}    },
    ${termination}
    nSims=${conditions_per_iteration},${balance}${constant}
)
<%def name="concentrations()">
% for spc in species_list:
    % if isinstance(spc["concentration"], (int, float)):
        '${spc["label"]}': ${spc["concentration"]},
    % endif
    % if isinstance(spc["concentration"], (tuple, list)):
        '${spc["label"]}': [${spc["concentration"][0]}, ${spc["concentration"][1]}],
    % endif
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
    % if isinstance(spc["concentration"], (int, float)):
        '${spc["label"]}': (${spc["concentration"]}, 'mol/cm^3'),
    % endif
    % if isinstance(spc["concentration"], (tuple, list)):
        '${spc["label"]}': [(${spc["concentration"][0]}, 'mol/cm^3'), (${spc["concentration"][1]}, 'mol/cm^3')],
    % endif
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
        if 'species_list' in reactor.keys():
            # This is relevant when a simulate adapter breaks ranged reactors down to individual conditions.
            species_list = reactor['species_list']
        else:
            # This is the base case when T3 generates an RMG input file for model generation.
            species_list = [{'label': spc['label'], 'concentration': spc['concentration']} for spc in species
                            if isinstance(spc['concentration'], (list, tuple))
                            or (isinstance(spc['concentration'], (float, int)) and spc['concentration'] > 0)
                            or spc['balance'] or not spc['reactive']]
            species_list.sort(key=lambda spc: spc['concentration'][0] if isinstance(spc['concentration'], (tuple, list))
                              else spc['concentration'], reverse=True)
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
                pressure = [(p, 'bar') for p in reactor['P']]
            else:
                raise ValueError(f"The reactor pressure must be a float or a list,\n"
                                 f"got {reactor['P']} which is a {type(reactor['P'])}.")
            balance = ''
            for spc in species:
                if spc['balance']:
                    balance = f"\n    balanceSpecies='{spc['label']}',"
                    break
            rmg_input += Template(gas_batch_constant_t_p_template).render(
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
    # Fall back to core_tolerance if tolerance_interrupt_simulation is None
    tis = model_input['tolerance_interrupt_simulation'] or model_input['core_tolerance']
    model['tolerance_interrupt_simulation'] = tis[iteration] \
        if len(tis) >= iteration + 1 else tis[-1]
    model_keys_to_skip = ['core_tolerance', 'tolerance_interrupt_simulation', 'atol', 'rtol', 'sens_atol', 'sens_rtol']
    args = ''
    for key, value in model_input.items():
        if key not in model_keys_to_skip and value is not None:
            args += f"\n    {to_camel_case(uv=key)}={value},"
    model['args'] = args
    rmg_input += Template(model_template).render(**model)

    # simulator
    if t3['sensitivity'] is not None and t3['sensitivity']['adapter'] == 'RMGConstantTP':
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


# Matches a 'method = ...' (or 'method=...') assignment statement, tolerant of both quote styles
# and of missing whitespace around '=': e.g. "    method = 'CSE',", '    method="CSE",',
# "method='CSE',". Anchored at the start of the (whitespace-stripped) line, so a line whose first
# non-whitespace character is '#' (a comment) never matches.
METHOD_LINE_RE = re.compile(r"^(?P<indent>[ \t]*)method\s*=\s*(?P<quote>['\"])(?P<value>.*?)(?P=quote)(?P<trailing>.*)$")

# A cheap, cheaper-than-the-full-rewrite-regex candidate check used by callers to decide whether a
# line is a 'method = ...' assignment at all (used to drive the per-file rewrite count). Kept
# separate from METHOD_LINE_RE so a malformed candidate (e.g. an assignment with an unterminated
# quote) still counts as "found" and fails loudly via rewrite_arkane_method_line's own ValueError,
# rather than being silently skipped as "not a method line".
METHOD_LINE_CANDIDATE_RE = re.compile(r'^[ \t]*method\s*=')


def rewrite_arkane_method_line(line: str, method: str) -> str:
    """
    Rewrite a P-dep network file's ``method = '...'`` line to use the resolved Arkane method name.

    Factored out of ``write_arkane_network_input_file`` so ``t3.pdep.hybrid`` can reuse the exact
    same string-scan logic rather than duplicating it.

    Tolerates both single- and double-quoted values (``method = 'CSE'`` and ``method = "CSE"``)
    and missing whitespace around ``=`` (``method='CSE'``), preserving whichever quote character
    the source line used. Does NOT attempt to distinguish a genuine top-level ``method = ...``
    assignment from one that happens to appear inside a triple-quoted string; that is a known,
    still-open limitation of this line-scan (non-AST) approach.

    Args:
        line (str): The source line, e.g. ``"    method = 'chemically-significant eigenvalues',\\n"``.
        method (str): 'CSE', 'MSC' or 'RS'.

    Raises:
        ValueError: If ``line`` does not have the expected ``method = <quoted string>`` shape.
            Never silently returns ``line`` unchanged: a caller that fails to rewrite a P-dep
            network's method line would go on to solve it with the SOURCE file's method while
            downstream cache metadata records the REQUESTED method, silently mismatching them.

    Returns:
        str: The rewritten line.
    """
    #     method = 'chemically-significant eigenvalues',
    newline = ''
    body = line
    if body.endswith('\r\n'):
        newline = '\r\n'
        body = body[:-2]
    elif body.endswith('\n'):
        newline = '\n'
        body = body[:-1]

    match = METHOD_LINE_RE.match(body)
    if match is None:
        raise ValueError(f"Could not rewrite a P-dep network file's 'method = ...' line to method {method!r}: the "
                         f"line does not have the expected quoted assignment shape (single- or double-quoted), "
                         f"got: {line!r}.")

    quote = match.group('quote')
    return (f"{match.group('indent')}method = {quote}{METHOD_MAP[method]}{quote}"
            f"{match.group('trailing')}{newline}")


def write_arkane_network_input_file(source_path: str,
                                    dest_path: str,
                                    method: str,
                                    sensitivity: bool = True,
                                    ) -> ArkaneNetworkWriteResult:
    """
    Rewrite an RMG P-dep network file into an Arkane network input file.

    Copies ``source_path`` to ``dest_path``, rewrites the ``method = '...'`` line via
    ``METHOD_MAP``, and (optionally) injects a ``sensitivity_conditions`` directive spanning the
    network's T/P extrema. This is the shared core behind both ``write_pdep_network_file`` (SA,
    ``sensitivity=True``) and the ``ArkaneMESolverAdapter`` (plain ME, ``sensitivity=False``).

    Args:
        source_path (str): The path to the original RMG network file to copy from.
        dest_path (str): The path to write the resulting Arkane input file to. The parent
                         directory is created if it does not already exist.
        method (str): 'CSE', 'MSC' or 'RS'.
        sensitivity (bool, optional): Whether to inject a ``sensitivity_conditions`` directive
                                      spanning the network's T/P extrema. Default: ``True``.

    Raises:
        ValueError: If T/P ranges could not be parsed from the file.

    Returns:
        ArkaneNetworkWriteResult: The current network's isomer labels, plus T-grid clamp
            provenance (see ``TGridClampRecord``).
    """
    dest_dir = os.path.dirname(dest_path)
    if dest_dir and not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)
    shutil.copyfile(src=source_path, dst=dest_path)

    with open(dest_path, 'r') as f:
        lines = f.readlines()
        # Computed once, from the exact text that was just copied to dest_path, BEFORE the T/P
        # extrema are line-scanned below: this is the ceiling every species the ceiling could be
        # computed for has valid NASA thermo up to, and it is what decides whether the pdep
        # block's own requested Tmax must be clamped. See network_thermo_t_max's docstring for
        # why None (no species contributed a determinable ceiling) means "clamp nothing" rather
        # than "clamp to zero", and for why a species whose thermo could not be read is skipped
        # rather than failing the whole file -- which is also why a non-empty `skipped` below
        # means this ceiling may be higher than the network's true one.
        ceiling = network_thermo_t_max(''.join(lines))
        thermo_t_max = ceiling.t_max
        if ceiling.skipped:
            logger.warning(
                f"Network '{source_path}': the NASA thermo ceiling used to clamp '{dest_path}' "
                f"({thermo_t_max} K) could not account for {len(ceiling.skipped)} species whose "
                f"thermo could not be read: {format_skipped_species(ceiling.skipped)}. The true "
                f"network-wide ceiling may therefore be lower than {thermo_t_max} K.")
        new_lines, isomer_labels = list(), list()
        t_min, t_max, t_count, p_min, p_max = None, None, None, None, None
        parse_tp, parse_isomers = False, (False, False)
        method_rewrite_count = 0
        # T-grid clamp provenance (see TGridClampRecord's docstring for why "clamped" is a
        # three-state design and why this is tracked at all): populated as the clamp decisions
        # below are actually made, not inferred after the fact from the final t_max/t_count.
        requested_t_max, clamped = None, False
        tlist_dropped, tlist_original_highest = False, None
        for line in lines:
            skip_line = False
            if 'pressureDependence(' in line:
                parse_tp = True
            if 'network(' in line:
                parse_isomers = (True, False)
            if parse_isomers[0] and 'isomers =' in line:
                parse_isomers = (True, True)
            if 'reactants =' in line:
                parse_isomers = (False, False)
            if parse_tp:
                if 'Tmin' in line and '(' in line:
                    #     Tmin = (300, 'K'),
                    t_min = line.split('(')[1].split(',')[0]
                elif 'Tmax' in line and '(' in line:
                    #     Tmax = (2200, 'K'),
                    t_max = line.split('(')[1].split(',')[0]
                    requested_t_max = float(t_max)
                    if thermo_t_max is not None and float(t_max) > thermo_t_max:
                        # The pdep grid asks for a temperature no species' NASA thermo is valid
                        # to; standalone Arkane refuses this outright ("No valid NASA polynomial
                        # at temperature ... K.") rather than tolerating the extrapolation the way
                        # RMG did at generation time. Clamp to the narrowest species ceiling, and
                        # rewrite THIS line in place (not just the sensitivity_conditions block
                        # below) since this is the line Arkane's pressureDependence job itself
                        # reads and fails on.
                        if t_min is None or thermo_t_max <= float(t_min):
                            raise ValueError(
                                f"Cannot write '{dest_path}': the pdep network's requested Tmax "
                                f"({t_max} K) exceeds the narrowest species thermo ceiling in the "
                                f"network ({thermo_t_max} K), and clamping to that ceiling would "
                                f"leave a Tmax ({thermo_t_max} K) at or below Tmin ({t_min} K), "
                                f"i.e. no valid temperature point to solve at all.")
                        logger.warning(
                            f"Network '{source_path}': requested pdep Tmax of {t_max} K exceeds "
                            f"the narrowest species thermo (NASA) ceiling in the file, "
                            f"{thermo_t_max} K; clamping the Arkane SA/ME grid's Tmax to "
                            f"{thermo_t_max} K in '{dest_path}' so standalone Arkane does not "
                            f"refuse the solve with 'No valid NASA polynomial at temperature ... K.'")
                        formatted_t_max = format_clamped_t_max(thermo_t_max)
                        paren_index = line.index('(')
                        comma_index = line.index(',', paren_index)
                        line = line[:paren_index + 1] + formatted_t_max + line[comma_index:]
                        t_max = formatted_t_max
                        clamped = True
                elif 'Tcount' in line and '=' in line:
                    #     Tcount = 8,
                    t_count = line.split('=', 1)[1].split(',')[0].strip()
                elif 'Tlist' in line and '(' in line:
                    #     Tlist = ([3200,2290.91,...,700],'K'),
                    # RMG's writer always precomputes an explicit Tlist alongside Tmin/Tmax/Tcount,
                    # and Arkane's own pdep.py only calls generate_T_list() (which would build the
                    # grid from Tmin/Tmax/Tcount) when self.Tlist is None. If Tlist is present, it
                    # wins outright and Tmin/Tmax are ignored -- so clamping the Tmax line above is
                    # a no-op unless this stale, out-of-range Tlist is also dropped so Arkane
                    # regenerates the grid over the clamped range instead of solving at the
                    # network's original (too-high) grid points.
                    if thermo_t_max is not None:
                        rhs = line.split('=', 1)[1].strip()
                        if rhs.endswith(','):
                            rhs = rhs[:-1]
                        try:
                            parsed = ast.literal_eval(rhs)
                        except (ValueError, SyntaxError) as e:
                            raise ValueError(
                                f"Cannot write '{dest_path}': a clamp is in play (the narrowest "
                                f"species thermo ceiling is {thermo_t_max} K), so the pdep "
                                f"network's Tlist line must be understood to decide whether to "
                                f"drop it, but it could not be parsed as a Python literal: {e}")
                        if (not isinstance(parsed, tuple) or len(parsed) != 2
                                or not isinstance(parsed[0], (list, tuple))
                                or not all(isinstance(v, (int, float)) and not isinstance(v, bool)
                                           for v in parsed[0])
                                or parsed[1] != 'K'):
                            raise ValueError(
                                f"Cannot write '{dest_path}': a clamp is in play (the narrowest "
                                f"species thermo ceiling is {thermo_t_max} K), so the pdep "
                                f"network's Tlist line must be understood to decide whether to "
                                f"drop it, but it is not a 2-element (sequence-of-numbers, unit) "
                                f"tuple with unit 'K'; got: {parsed!r}.")
                        entries = parsed[0]
                        if any(entry > thermo_t_max for entry in entries):
                            missing = [name for name, val in (('Tmin', t_min), ('Tmax', t_max),
                                                              ('Tcount', t_count)) if val is None]
                            if missing:
                                raise ValueError(
                                    f"Cannot write '{dest_path}': the pdep network's Tlist contains "
                                    f"an entry above the narrowest species thermo ceiling "
                                    f"({thermo_t_max} K) and must be dropped so Arkane regenerates "
                                    f"the T grid itself, but {', '.join(missing)} could not be "
                                    f"parsed from the pressureDependence(...) block to regenerate "
                                    f"from.")
                            highest = max(entries)
                            logger.warning(
                                f"Network '{source_path}': dropping the pressureDependence(...) "
                                f"block's explicit Tlist line when writing '{dest_path}' -- its "
                                f"highest entry ({highest} K) exceeds the narrowest species thermo "
                                f"(NASA) ceiling in the file, {thermo_t_max} K; Arkane will "
                                f"regenerate the T grid itself from the clamped Tmin ({t_min} K) / "
                                f"Tmax ({t_max} K) / Tcount ({t_count}) instead of using this "
                                f"network's original (too-high) grid.")
                            skip_line = True
                            tlist_dropped = True
                            tlist_original_highest = float(highest)
                elif 'Pmin' in line and '(' in line:
                    #     Pmin = (0.01, 'bar'),
                    p_min = line.split('(')[1].split(',')[0]
                elif 'Pmax' in line and '(' in line:
                    #     Pmax = (100, 'bar'),
                    p_max = line.split('(')[1].split(',')[0]
            if all(parse_isomers) and "'," in line and "'" in line:
                #         'C=O(26)',
                parts = line.split("'")
                if len(parts) >= 2:
                    isomer_labels.append(parts[1])
            if METHOD_LINE_CANDIDATE_RE.match(line):
                new_lines.append(rewrite_arkane_method_line(line=line, method=method))
                method_rewrite_count += 1
            elif 'rmgmode' in line:
                new_lines.append(line)
                if any(param is None for param in [t_min, t_max, p_min, p_max]):
                    raise ValueError(f'Could not parse all T/P parameters, got:\n'
                                     f'T min = {t_min}, T max = {t_max}, P min = {p_min}, P max = {p_max}.')
                if sensitivity:
                    sa_conditions = f"""    sensitivity_conditions = [[({t_min}, 'K'), ({p_min}, 'bar')],
                              [({t_max}, 'K'), ({p_min}, 'bar')],
                              [({t_min}, 'K'), ({p_max}, 'bar')],
                              [({t_max}, 'K'), ({p_max}, 'bar')]],"""
                    new_lines.append(sa_conditions)
            elif not skip_line:
                new_lines.append(line)

    if method_rewrite_count != 1:
        raise ValueError(f"Expected exactly one 'method = ...' line to rewrite in '{source_path}', found "
                         f"{method_rewrite_count}. Refusing to write '{dest_path}' with an unresolved (zero) or "
                         f"ambiguous (more than one) method line, rather than silently solving with the source "
                         f"file's original method.")

    with open(dest_path, 'w') as f:
        f.writelines(new_lines)

    t_grid_clamp = TGridClampRecord(
        clamped=clamped,
        requested_t_max=requested_t_max,
        thermo_ceiling=thermo_t_max,
        written_t_max=float(t_max) if t_max is not None else None,
        tlist_dropped=tlist_dropped,
        tlist_original_highest=tlist_original_highest,
        skipped_species=tuple(ceiling.skipped),
    )
    return ArkaneNetworkWriteResult(isomer_labels=tuple(isomer_labels), t_grid_clamp=t_grid_clamp)


def write_pdep_network_file(network_name: str,
                            method: str,
                            pdep_sa_path: str,
                            rmg_pdep_path: str,
                            ) -> ArkaneNetworkWriteResult:
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
        ArkaneNetworkWriteResult: The current network's isomer labels, plus T-grid clamp
            provenance (see ``TGridClampRecord``).
    """
    sa_pdep_path = os.path.join(pdep_sa_path, network_name, method)
    input_file_path = os.path.join(sa_pdep_path, 'input.py')
    source_path = os.path.join(rmg_pdep_path, network_name + '.py')
    return write_arkane_network_input_file(source_path=source_path,
                                           dest_path=input_file_path,
                                           method=method,
                                           sensitivity=True)


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


def get_species_obj_from_a_species_dict(species_dict: dict,
                                        raise_error: bool = False,
                                        ) -> T3Species | None:
    """
    Get a T3Species instance that corresponds to a species specified under the rmg.species
    section of the T3 input file (a species dictionary).

    Args:
        species_dict (dict): The species dictionary to process.
        raise_error (bool, optional): Whether to raise an error if a T3Species instance cannot be generated.
                                      Default: ``False``.

    Raises:
        ValueError: If the species dictionary does not have a specified structure (if ``raise_error`` is ``True``).

    Returns:
        T3Species: The corresponding T3Species instance, or ``None`` if it cannot be generated.
    """
    species, errored = None, False
    if species_dict['adjlist'] is not None:
        species = T3Species(label=species_dict['label'], adjlist=species_dict['adjlist'])
    elif species_dict['smiles'] is not None:
        species = T3Species(label=species_dict['label'], smiles=species_dict['smiles'])
    elif species_dict['inchi'] is not None:
        species = T3Species(label=species_dict['label'], inchi=species_dict['inchi'])
    elif species_dict['xyz'] is not None:
        xyz_entries = species_dict['xyz'] if isinstance(species_dict['xyz'], list) else [species_dict['xyz']]
        for xyz in xyz_entries:
            mol = perceive_molecule_from_xyz(xyz=xyz)
            if mol is not None:
                species = T3Species(label=species_dict['label'], mol=mol)
                break
        else:
            errored = True
    else:
        errored = True
    if errored and raise_error:
        raise ValueError(f"The species corresponding to {species_dict['label']} does not have a specified structure.")
    return species
