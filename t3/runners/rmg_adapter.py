
import os
import math

from mako.template import Template
from t3.utils.writer import to_camel_case
from t3.common import get_rmg_species_from_a_species_dict
from t3.utils.generator import generate_radicals
from t3.imports import local_t3_path, settings, submit_scripts
from t3.utils.ssh import SSHClient
from t3.logger import Logger


METHOD_MAP = {'CSE': 'chemically-significant eigenvalues',
              'RS': 'reservoir state',
              'MSC': 'modified strong collision',
              }

RMG_EXECUTION_TYPE = settings['execution_type']['rmg']
submit_filenames = settings['submit_filenames']
rmg_memory = settings['rmg_initial_memory']
if RMG_EXECUTION_TYPE == 'queue':
    SERVER = list(settings['servers'].keys())[0]
    CPUS = settings['servers'][SERVER]['cpus']
    MEMORY = settings['servers'][SERVER]['memory']
    CLUSTER_SOFT = settings['servers'][SERVER]['cluster_soft']
elif RMG_EXECUTION_TYPE == 'local':
    SERVER = 'local'
elif RMG_EXECUTION_TYPE == 'incore':
    SERVER = 'incore'
else:
    raise ValueError(f'RMG execution type {RMG_EXECUTION_TYPE} is not supported.')



class RMGAdapter(object):
    
    def __init__(self, 
                 rmg: dict,
                 t3: dict,
                 iteration: int,
                 paths: dict,
                 walltime: str="00:00:00",
                 cpus: int=1,
                 memory: str="8G", # TODO: make this a parameter
                 max_iterations: int=1,
                 verbose: bool=False,
                 t3_project_name: str=None,
                 rmg_execution_type: str='incore',
                 restart_rmg: bool=False,
                 server: str=None,                 
                 ):
        self.rmg = rmg
        self.t3 = t3
        self.iteration = iteration
        self.paths = paths
        self.walltime = walltime
        self.rmg_path = self.paths['RMG']
        self.rmg_input_file_path = self.paths['RMG input']
        self.max_cpus = CPUS if CPUS else cpus
        self.max_memory = MEMORY if MEMORY else memory
        self.t3_project_name = t3_project_name
        self.max_iterations = max_iterations
        self.rmg_execution_type = RMG_EXECUTION_TYPE or rmg_execution_type
        if self.rmg_execution_type == 'queue':
            self.max_job_time = settings['servers'][SERVER]['max_job_time']
        
        
        
        
        if not os.path.isdir(local_t3_path):
            os.makedirs(local_t3_path)
        
        self.files_to_upload = list()
        self.folder_to_download = None
        
    def run_rmg(self):
        """
        Run RMG
        """
        self.set_cpu_and_mem()
        self.set_file_paths()
        self.set_files()

        if self.rmg_execution_type == 'incore':
            self.execute_incore()
        elif self.rmg_execution_type == 'queue':
            self.execute_queue()
        else:
            raise ValueError(f'RMG execution type {self.rmg_execution_type} is not supported.')
        
    def write_rmg_input_file(self):
        """
        Write an RMG input file to the given file path.
        Will create the directory if needed.
        """
        rmg= self.rmg.copy()
        rmg_input = ''
        iteration = self.iteration - 1 # iteration is 1-indexed, convert to 0-indexed for list indexing
        
        # Database
        database = rmg['database']
        # The following args types could be either str or list, detect str and format accordingly
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
                raise ValueError(f"A species must have either an adjlist, smiles, or inchi. Species {spc['label']} has none of these.")
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
        
        # Solvent
        solvent_template = """solvation(solvent='${solvent}')

"""
        solvent = ''
        for spc in species:
            # The schema assures that there is only one species define as a solvent
            # TODO: Assue that the requested solvent actually exists in the RMG database
            if spc['solvent']:
                solvent = spc['label']
                break
        
        if solvent:
            rmg_input += Template(solvent_template).render(solvent=solvent)
        
        # Model
        model_input = rmg['model']
        model_template = """model(
    toleranceMoveToCore=${tol_move_to_core},
    toleranceInterruptSimulation=${tolerance_interrupt_simulation},${args}
)
"""
        model = dict()
        model['tol_move_to_core'] = model_input['core_tolerance'][iteration]\
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
        
        # Simulator
        if self.t3['sensitivity'] is not None and self.t3['sensitivity']['adapter'] == 'RMGConstantTP':
            simulator_template = """\nsimulator(atol=${atol}, rtol=${rtol}, sens_atol=${sens_atol}, sens_rtol=${sens_rtol})\n"""
            rmg_input += Template(simulator_template).render(atol=model_input['atol'],
                                                            rtol=model_input['rtol'],
                                                            sens_atol=self.t3['sensitivity']['atol'],
                                                            sens_rtol=self.t3['sensitivity']['rtol']
                                                            )
        else:
            simulator_template = """\nsimulator(atol=${atol}, rtol=${rtol})\n"""
            rmg_input += Template(simulator_template).render(atol=model_input['atol'],
                                                            rtol=model_input['rtol'],
                                                            )
        
        # PressureDependence
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
            options['walltime'] = self.walltime
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

        if not os.path.isdir(os.path.dirname(self.rmg_input_file_path)):
            os.makedirs(os.path.dirname(self.rmg_input_file_path))
        with open(self.rmg_input_file_path, 'w') as f:
            f.writelines(rmg_input)

    def set_cpu_and_mem(self):
        """
        Set cpu and memory based on ESS and cluster software.
        This is not an abstract method and should not be overwritten.
        """
        if self.max_memory < rmg_memory:
            Logger.warning(f'The maximum memory of the server is {self.max_memory} GB, which is less than the defined '
                           f'RMG memory {rmg_memory} GB set by the user. The RMG memory will be set to {self.max_memory}')
            total_submit_script_memory = self.max_memory
        else:
            total_submit_script_memory = rmg_memory
        # Determine amount of memory in submit script based on cluster job scheduling system.
        cluster_software = CLUSTER_SOFT.lower() if SERVER is not None else None
        if cluster_software in ['oge', 'sge', 'htcondor']:
            # In SGE, "-l h_vmem=5000M" specifies the memory for all cores to be 5000 MB.
            self.submit_script_memory = math.ceil(total_submit_script_memory)  # in MB
        if cluster_software in ['pbs']:
            # In PBS, "#PBS -l select=1:ncpus=8:mem=12000000" specifies the memory for all cores to be 12 MB.
            self.submit_script_memory = math.ceil(total_submit_script_memory) * 1E6  # in Bytes
        elif cluster_software in ['slurm']:
            # In Slurm, "#SBATCH --mem-per-cpu=2000" specifies the memory **per cpu/thread** to be 2000 MB.
            self.submit_script_memory = math.ceil(total_submit_script_memory / self.max_cpus)  # in MB
        self.set_input_file_memory()


    def set_files(self) -> None:
        """
        Set files to be uploaded and downloaded. Writes the files if needed.
        Modifies the self.files_to_upload and self.files_to_download attributes.

        self.files_to_download is a list of remote paths.

        self.files_to_upload is a list of dictionaries, each with the following keys:
        ``'name'``, ``'source'``, ``'make_x'``, ``'local'``, and ``'remote'``.
        If ``'source'`` = ``'path'``, then the value in ``'local'`` is treated as a file path.
        Else if ``'source'`` = ``'input_files'``, then the value in ``'local'`` will be taken
        from the respective entry in inputs.py
        If ``'make_x'`` is ``True``, the file will be made executable.
        """
        # 1. ** Upload **
        # 1.1. Submit script
        if self.rmg_execution_type != 'incore':
            # We need a submit script (submitted local or SSH)
            self.write_submit_script()
            self.files_to_upload.append(self.get_file_property_dictionary(
                file_name='submit.sh'))
        # 1.2. RMG input file
        self.write_rmg_input_file()
        self.files_to_upload.append(self.get_file_property_dictionary(file_name='input.py'))
        
        # 2. ** Download **
        # 2.1. RMG output folder
        self.folder_to_download = self.remote_path
        
    
    def set_additional_file_paths(self) -> None:
        """
        Set additional file paths to be uploaded and downloaded.
        Modifies the self.additional_file_paths attribute.
        """
        pass
    
    def set_input_file_memory(self) -> None:
        """
        Set the self.input_file_memory attribute.
        """
        pass
    
    def execute_incore(self) -> None:
        """
        Execute the job incore.
        """
        pass
    
    def execute_queue(self) -> None:
        """
        Execute the job in the queue.
        """
        if SERVER != 'local':
            with SSHClient(SERVER) as ssh:
                self.job_status, self.job_id = ssh.submit_job(remote_path=self.remote_path)
    
    def execute_local(self) -> None:
        """
        Execute the job.
        """
        if SERVER == 'local':
            self.job_status, self.job_id = self.submit_job()

    def write_submit_script(self) -> None:
        """
        Write the submit script.
        """
        if SERVER is None:
            return
        if self.max_job_time < self.walltime:
            self.walltime = self.max_job_time
        architecture = ''
        
        try:
            submit_script = submit_scripts['rmg'][settings['servers'][SERVER]['cluster_soft']].format(
                                            name=f'{self.t3_project_name}_RMG_{self.iteration}',
                                            max_job_time=self.max_job_time,
                                            cpus=self.max_cpus,
                                            memory=self.submit_script_memory,
                                            max_iterations=self.max_iterations,
                                            )
        except KeyError as e:
            raise KeyError(f'Invalid key in submit script: {e}')
    
        with open(os.path.join(self.rmg_path, 'submit.sh'), 'w') as f:
            f.write(submit_script)
            
    def set_file_paths(self) -> None:
        """
        Set local and remote file paths.
        """
        
        self.local_iteration_path = self.paths['iteration']
        self.local_rmg_path = self.paths['RMG']
        
        
        if SERVER != 'incore':
            path = settings['servers'][SERVER].get('path','').lower()
            path = os.path.join(path, settings['servers'][SERVER]['un']) if path else ''
            self.remote_path = os.path.join(path, 'runs', 'T3_Projects', self.t3_project_name, f"iteration_{self.iteration}", 'RMG')
        
        # Get additional file paths - but I think we copy the whole folder of RMG from Remote to Local
        self.set_additional_file_paths()
        
    def get_file_property_dictionary(self,
                                     file_name: str,
                                     local: str='',
                                     remote: str='',
                                     source: str='path',
                                     make_x: bool=False) -> dict:
        """
        Get a dictionary that represents a file to be uploaded or downloaded to/from a server via SSH.

        Args:
            file_name (str): The file name.
            local (str, optional): The full local path.
            remote (str, optional): The full remote path.
            source (str, optional): Either ``'path'`` to treat the ``'local'`` attribute as a file path,
                                    or ``'input_files'`` to take the respective entry from inputs.py.
            make_x (bool, optional): Whether to make the file executable, default: ``False``.

        Returns:
            dict: A file representation.
        """
        if not file_name:
            raise ValueError('file_name must be specified')
        if source not in ['path', 'input_files']:
            raise ValueError(f'source must be either "path" or "input_files", got {source}')
        local = local or os.path.join(self.local_rmg_path, file_name)
        remote = remote or os.path.join(self.remote_path, file_name) if self.remote_path else None
        return{'file_name': file_name,
               'local': local,
               'remote': remote,
               'source': source,
               'make_x': make_x}
    