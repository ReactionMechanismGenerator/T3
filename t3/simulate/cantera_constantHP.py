"""
Cantera Simulator Adapter module
Used to run mechanism analysis with Cantera as an ideal gas in a batch reactor at constant H-P
"""

import cantera as ct
import numpy as np
from typing import List, Optional, Type

from rmgpy.tools.canteramodel import generate_cantera_conditions
from rmgpy.tools.data import GenericData

from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter


class CanteraConstantHP(SimulateAdapter):
    """
    CanteraConstantHP is an adapter for the abstract class SimulateAdapter that simulates ideal gases
    in a batch reactor at constant pressure and constant enthalpy.

    Args:
        t3 (dict): The T3.t3 attribute, which is a dictionary containing the t3 block from the input yaml or API.
        rmg (dict): The T3.rmg attribute, which is a dictionary containing the rmg block from the input yaml or API.
        paths (dict): The T3.paths attribute, which is a dictionary containing relevant paths.
        logger (Logger): Instance of T3's Logger class.
        atol (float): The absolute tolerance used when integrating.
        rtol (float): The relative tolerance used when integrating.
        observable_list (Optional[list]): Species used for SA. Entries are species labels as strings. Example: ['OH']
        sa_atol (float, optional): The absolute tolerance used when performing sensitivity analysis.
        sa_atol (float, optional): The relative tolerance used when performing sensitivity analysis.
        global_observables (Optional[List[str]]): List of global observables ['IgD', 'ESR', 'SL'] used by Cantera adapters.

    Attributes:
        all_data (list): List containing the following RMG GenericData objects grouped as a tuple:
                         time, data_list, reaction_sensitivity_data, thermodynamic_sensitivity_data
        atol (float): The absolute tolerance used when integrating.
        cantera_reactor_type (str): String specifying the type of Cantera reactor to use.
        cantera_simulation (ct.ReactorNet): Cantera reactor net object.
        conditions (list): List whose entries are reaction conditions for simulation.
        global_observables (List[str]): List of global observables ['IgD', 'ESR', 'SL'] used by Cantera adapters.
        inert_list (list): List of possible inert species in the model
        inert_index_list (list): List of indices corresponding to the inert species present in the model.
        initialconds (dict): Key is the Cantera species. Value is the initial mol fraction.
        logger (Logger): Instance of T3's Logger class.
        model (ct.Solution): Cantera solution object for the mechanism.
        num_ct_reactions (int): Number of reactions in the model.
        num_ct_species (int): Number of species in the model.
        observable_list (list): Species used for SA. Entries are species labels as strings. Example: ['OH']
        paths (dict): The T3.paths attribute, which is a dictionary containing relevant paths.
        rmg (dict): The T3.rmg attribute, which is a dictionary containing the rmg block from the input yaml or API.
        rtol (float): The relative tolerance used when integrating.
        rxn_identifier_lookup (dict): Keys are reactions (str). Values are index in the model.
        sa_atol (float): The absolute tolerance used when performing sensitivity analysis.
        sa_rtol (float): The relative tolerance used when performing sensitivity analysis.
        sensitive_species (list): List of sensitive species. Entries are strings that include the RMG index.
        spc_identifier_lookup (dict): Keys are species (str). Values are index in the model.
        t3 (dict): The T3.t3 attribute, which is a dictionary containing the t3 block from the input yaml or API.
    """

    def __init__(self,
                 t3: dict,
                 rmg: dict,
                 paths: dict,
                 logger: Type[Logger],
                 atol: float = 1e-16,
                 rtol: float = 1e-8,
                 observable_list: Optional[list] = None,
                 sa_atol: float = 1e-6,
                 sa_rtol: float = 1e-4,
                 global_observables: Optional[Type[List[str]]] = None
                 ):

        self.t3 = t3
        self.rmg = rmg
        self.paths = paths
        self.logger = logger
        self.atol = atol
        self.rtol = rtol
        self.observable_list = observable_list or list()
        self.sa_atol = sa_atol
        self.sa_rtol = sa_rtol
        self.global_observables = global_observables

        # initialize other attributes
        self.sensitive_species = list()
        self.initialconds = dict()
        self.all_data = list()
        # this adapter is for constant P batch simulations
        self.cantera_reactor_type = 'IdealGasConstPressureReactor'

        self.model = None
        self.cantera_simulation = None
        self.conditions = list()
        self.inert_list = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'N2']
        self.inert_index_list = list()

        self.num_ct_reactions = None
        self.num_ct_species = None

        self.set_up()

        # create list of indices of the inerts
        for i, species in enumerate(self.model.species()):
            if species.name in self.inert_list:
                self.inert_index_list.append(i)

        self.spc_identifier_lookup = {}
        for i, spc in enumerate(self.model.species()):
            self.spc_identifier_lookup[spc.name] = i

        self.rxn_identifier_lookup = {}
        for i, rxn in enumerate(self.model.reactions()):
            self.rxn_identifier_lookup[rxn.equation] = i

    def set_up(self):
        """
        Read in the Cantera input file and set up the empty attributes initialized in the init method.
        """
        # read in the cantera .cti file
        self.model = ct.Solution(infile=self.paths['cantera annotated'])
        self.num_ct_reactions = len(self.model.reactions())
        self.num_ct_species = len(self.model.species())

        # create list of indices of the inerts
        for i, species in enumerate(self.model.species()):
            if species.name in self.inert_list:
                self.inert_index_list.append(i)

        self.spc_identifier_lookup = {}
        for i, spc in enumerate(self.model.species()):
            self.spc_identifier_lookup[spc.name] = i

        self.rxn_identifier_lookup = {}
        for i, rxn in enumerate(self.model.reactions()):
            self.rxn_identifier_lookup[rxn.equation] = i

        self.species_names_without_indices = [self.model.species()[i].name.split('(')[0] for i in
                                              range(self.num_ct_species)]
        # set initial conditions and find any species for SA
        for input_species in self.rmg['species']:
            # find index of this species in the list of Cantera species
            idx = self.species_names_without_indices.index(input_species['label'])
            self.initialconds.update({self.model.species(idx).name: input_species['concentration']})

            if self.species_names_without_indices[idx] in self.observable_list:
                self.sensitive_species.append(self.model.species(idx).name)

        reactor_type_list = [self.cantera_reactor_type]
        mol_frac_list = [self.initialconds]
        Tlist = ([self.rmg['reactors'][0]['T']], 'K')
        Plist = ([self.rmg['reactors'][0]['P']], 'bar')
        rxn_time, units = self.rmg['reactors'][0]['termination_time']
        reaction_time_list = ([rxn_time], units)  # tuple giving the ([list of reaction times], units)

        self.generate_conditions(reactor_type_list, reaction_time_list, mol_frac_list, Tlist, Plist)

    def generate_conditions(self,
                            reactor_type_list: Type[List[tuple]],
                            reaction_time_list: Type[List[tuple]],
                            mol_frac_list: Type[List[dict]],
                            T0_list: Type[List[tuple]] = None,
                            P0_list: Type[List[tuple]] = None,
                            V0_list: Type[List[tuple]] = None,
                            ):
        """
        Saves all the reaction conditions.

        Args:
            reactor_type_list (list): A list of strings specifying the type of supported Cantera reactor:
                                      IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
                                      IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures
                                      IdealGasConstPressureTemperatureReactor: A homogenous, constant pressure and constant temperature, zero-dimensional reactor
                                                                              for ideal gas mixtures (the same as RMG's SimpleReactor)
            reaction_time_list (tuple): A tuple object giving the ([list of reaction times], units)
            mol_frac_list (list): A list of molfrac dictionaries with species object keys
                               and mole fraction values

            To specify the system for an ideal gas, 2 of the following 3 parameters must be defined:
                T0_list (tuple): A tuple giving the ([list of initial temperatures], units)
                P0_list (tuple): A tuple giving the ([list of initial pressures], units)
                V0_list (tuple): A tuple giving the ([list of initial specific volumes], units)
        """

        self.conditions = generate_cantera_conditions(reactor_type_list,
                                                      reaction_time_list,
                                                      mol_frac_list,
                                                      T0_list,
                                                      P0_list,
                                                      V0_list,
                                                      )

    def reinitialize_simulation(self,
                                T0=None,
                                P0=None,
                                X0=None,
                                V0=None,
                                ):
        """
        Re-initializes the cantera solution object (self.model) and cantera reactor object (self.cantera_reactor).
        This method is called at the start of other methods in this class.

        Args:
            T0 (float): Initial temperature in Kelvin.
            P0 (float): Initial pressure in Pascals.
            X0 (dict): Initial mole fractions.
            V0 (float): Initial volume in m3.
        """
        # assign initial conditions
        if V0 is None:
            self.model.TPX = T0, P0, X0
        elif P0 is None:
            self.model.TDX = T0, 1 / V0, X0

        self.cantera_reactor = ct.IdealGasConstPressureReactor(contents=self.model)
        # Run this individual condition as a simulation
        self.cantera_simulation = ct.ReactorNet([self.cantera_reactor])

        # set simulation tolerance
        self.cantera_simulation.atol = self.atol
        self.cantera_simulation.rtol = self.rtol

        if self.sensitive_species:
            if ct.__version__ == '2.2.1':
                self.logger.warning('Cantera version 2.2.1 may not support sensitivity analysis unless SUNDIALS '
                                    'was used during compilation.')
                self.logger.warning('Upgrade to newer of Cantera in anaconda using the command '
                                    '"conda update -c rmg cantera"')
            # Add all the reactions as part of SA
            for i in range(self.num_ct_reactions):
                self.cantera_reactor.add_sensitivity_reaction(i)
            # add all species enthalpies as part of SA
            for i in range(self.num_ct_species):
                self.cantera_reactor.add_sensitivity_species_enthalpy(i)
            # Set the tolerances for the sensitivity coefficients
            self.cantera_simulation.rtol_sensitivity = self.sa_atol
            self.cantera_simulation.atol_sensitivity = self.sa_rtol

    def simulate(self):
        """
        Simulate the mechanism and store all results to the all_data attribute.
        """

        if self.sensitive_species:
            self.logger.info('Running a simulation with SA using CanteraConstantHP...')
        else:
            self.logger.info('Running a simulation using CanteraConstantHP...')

        species_names_list = [species.name for species in self.model.species()]
        self.all_data = list()

        for condition in self.conditions:
            # Set Cantera simulation conditions
            T0 = condition.T0.value_si
            try:
                V0 = self.conditions[0].V0.value_si
                P0 = None
            except AttributeError as e:
                P0 = condition.P0.value_si
                V0 = None
            self.reinitialize_simulation(T0=T0,
                                         P0=P0,
                                         X0=condition.mol_frac,
                                         V0=V0,
                                         )

            # Initialize the variables to be saved
            times = []
            temperature = []
            pressure = []
            species_data = []
            kinetic_sensitivity_data = []
            thermo_sensitivity_data = []

            # Begin integration
            while self.cantera_simulation.time < condition.reaction_time.value_si:

                # Advance the state of the reactor network in time from the current time to time t [s], taking as many integrator timesteps as necessary.
                self.cantera_simulation.step()
                times.append(self.cantera_simulation.time)
                temperature.append(self.cantera_reactor.T)
                pressure.append(self.cantera_reactor.thermo.P)
                species_data.append(self.cantera_reactor.thermo[species_names_list].X)

                if self.sensitive_species:
                    # Cantera returns mass-based sensitivities rather than molar concentration or mole fraction based sensitivities.
                    # The equation for converting between them is:
                    #
                    # d ln xi = d ln wi - sum_(species i) (dln wi) (xi)
                    #
                    # where xi is the mole fraction of species i and wi is the mass fraction of species i

                    mass_frac_sensitivity_array = self.cantera_simulation.sensitivities()
                    if condition.reactor_type == 'IdealGasReactor':
                        # Row 0: mass, Row 1: volume, Row 2: internal energy or temperature, Row 3+: mass fractions of species
                        mass_frac_sensitivity_array = mass_frac_sensitivity_array[3:, :]
                    elif condition.reactor_type == 'IdealGasConstPressureReactor' or condition.reactor_type == 'IdealGasConstPressureTemperatureReactor':
                        # Row 0: mass, Row 1: enthalpy or temperature, Row 2+: mass fractions of the species
                        mass_frac_sensitivity_array = mass_frac_sensitivity_array[2:, :]
                    else:
                        raise Exception('Other types of reactor conditions are currently not supported')

                    for i in range(len(mass_frac_sensitivity_array)):
                        mass_frac_sensitivity_array[i] *= species_data[-1][i]

                    # extract kinetics SA
                    kinetics_mass_frac_sa = mass_frac_sensitivity_array[:, 0:self.num_ct_reactions]
                    sensitivity_array = np.zeros(len(self.sensitive_species) * len(self.model.reactions()))
                    for index, species in enumerate(self.sensitive_species):
                        for j in range(self.num_ct_reactions):
                            sensitivity_array[self.num_ct_reactions * index + j] = self.cantera_simulation.sensitivity(
                                species, j)

                            for i in range(len(kinetics_mass_frac_sa)):
                                if i not in self.inert_index_list:
                                    # massFracSensitivity for inerts are returned as 0.0 in Cantera, so we do not include them here
                                    sensitivity_array[self.num_ct_reactions * index + j] -= kinetics_mass_frac_sa[i][j]
                    kinetic_sensitivity_data.append(sensitivity_array)

                    # extract thermo SA
                    thermo_mass_frac_sa = mass_frac_sensitivity_array[:, self.num_ct_reactions:]
                    sensitivity_array = np.zeros(len(self.sensitive_species) * self.num_ct_species)
                    for index, species in enumerate(self.sensitive_species):
                        for j in range(self.num_ct_species):
                            sensitivity_array[self.num_ct_species * index + j] = self.cantera_simulation.sensitivity(
                                species, j + self.num_ct_reactions)

                            for i in range(len(mass_frac_sensitivity_array)):
                                if i not in self.inert_index_list:
                                    # massFracSensitivity for inerts are returned as 0.0 in Cantera, so we must not include them here
                                    sensitivity_array[self.num_ct_species * index + j] -= thermo_mass_frac_sa[i][j]
                    thermo_sensitivity_data.append(sensitivity_array)

            # Convert species_data and sensitivity data to numpy arrays
            species_data = np.array(species_data)
            kinetic_sensitivity_data = np.array(kinetic_sensitivity_data)
            thermo_sensitivity_data = np.array(thermo_sensitivity_data)

            # Resave data into generic data objects
            time = GenericData(label='Time',
                               data=times,
                               units='s')
            temperature = GenericData(label='Temperature',
                                      data=temperature,
                                      units='K')
            pressure = GenericData(label='Pressure',
                                   data=pressure,
                                   units='Pa')
            condition_data = []
            condition_data.append(temperature)
            condition_data.append(pressure)

            for index, species in enumerate(self.model.species()):
                # Create generic data object that saves the species object into the species object.  To allow easier manipulate later.
                species_generic_data = GenericData(label=species.name,
                                                   species=species,
                                                   data=species_data[:, index],
                                                   index=index
                                                   )
                condition_data.append(species_generic_data)

            # save kinetic data as generic data object
            reaction_sensitivity_data = []
            for index, species in enumerate(self.sensitive_species):
                for j in range(self.num_ct_reactions):
                    reaction_sensitivity_generic_data = GenericData(
                        label='dln[{0}]/dln[k{1}]: {2}'.format(species, j + 1, self.model.reactions()[j]),
                        species=species,
                        reaction=self.model.reactions()[j],
                        data=kinetic_sensitivity_data[:, self.num_ct_reactions * index + j],
                        index=j + 1,
                    )
                    reaction_sensitivity_data.append(reaction_sensitivity_generic_data)

            # save thermo data as generic data object
            thermodynamic_sensitivity_data = []
            for index, species in enumerate(self.sensitive_species):
                for j in range(self.num_ct_species):
                    thermo_sensitivity_generic_data = GenericData(
                        label='dln[{0}]/dH[{1}]'.format(species, self.model.species()[j].name),
                        species=species,
                        data=thermo_sensitivity_data[:, self.num_ct_species * index + j],
                        index=j + 1,
                    )
                    thermodynamic_sensitivity_data.append(thermo_sensitivity_generic_data)

            self.all_data.append((time, condition_data, reaction_sensitivity_data, thermodynamic_sensitivity_data))

    def get_sa_coefficients(self):
        """
        Obtain the SA coefficients.

        Returns:
             sa_dict (dict): a SA dictionary, whose structure is given in the docstring for T3/t3/main.py
        """
        sa_dict = {'kinetics': dict(), 'thermo': dict(), 'time': list()}

        for condition_data in self.all_data:
            time, data_list, reaction_sensitivity_data, thermodynamic_sensitivity_data = condition_data
            sa_dict['time'] = time.data

            # extract kinetic SA
            for rxn in reaction_sensitivity_data:
                # for kinetics, get `ethane(1)` from `dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)`
                observable_label = rxn.label.split('[')[1].split(']')[0]
                if observable_label not in sa_dict['kinetics']:
                    sa_dict['kinetics'][observable_label] = dict()
                # for kinetics, get k8 from `dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)` then only extract 8
                parameter = rxn.label.split('[')[2].split(']')[0]
                parameter = int(parameter[1:])
                sa_dict['kinetics'][observable_label][parameter] = rxn.data

            # extract thermo SA
            for spc in thermodynamic_sensitivity_data:
                # for thermo get 'C2H4(8)' from `dln[ethane(1)]/dH[C2H4(8)]`
                observable_label = spc.label.split('[')[1].split(']')[0]
                if observable_label not in sa_dict['thermo']:
                    sa_dict['thermo'][observable_label] = dict()
                # for thermo get 'C2H4(8)' from `dln[ethane(1)]/dH[C2H4(8)]`
                parameter = spc.label.split('[')[2].split(']')[0]
                sa_dict['thermo'][observable_label][parameter] = spc.data

        return sa_dict

    def get_idt_by_T(self):
        """
        Finds the ignition point by approximating dT/dt as a first order forward difference
        and then finds the point of maximum slope.

        Returns:
            idt_dict (dict): Dictionary whose keys are 'idt' and 'idt_index' and whose values are lists of
                             the ignition delay time in seconds and index at which idt occurs respectively.
        """
        idt_dict = {'idt': list(),
                    'idt_index': list(),
                    }
        for i, condition_data in enumerate(self.all_data):
            time, data_list, reaction_sensitivity_data, thermodynamic_sensitivity_data = condition_data
            T_data = data_list[0]

            dTdt = np.diff(T_data.data) / np.diff(time.data)
            idt_dict['idt_index'].append(int(np.argmax(dTdt)))
            idt_dict['idt'].append(time.data[idt_dict['idt_index'][i]])

        return idt_dict

    def find_equilibrium(self,
                         constrained_state_vars: str,
                         ):
        """
        Args:
            constrained_state_vars (str): One of the following options supported by Cantera:
                                           'TP','TV','HP','SP','SV','UV'

        Returns:
              equilibrium_dictionaries (list[dict]): List whose entries are the mole fraction dictionaries of the
                                                     equilibrated model.
        """
        equilibrium_dictionaries = list()
        for condition in self.conditions:
            # Set Cantera simulation conditions
            T0 = condition.T0.value_si
            try:
                V0 = self.conditions[0].V0.value_si
                P0 = None
            except AttributeError as e:
                P0 = condition.P0.value_si
                V0 = None

            self.reinitialize_simulation(T0=T0,
                                         P0=P0,
                                         X0=condition.mol_frac,
                                         V0=V0,
                                         )

            self.model.equilibrate(constrained_state_vars)
            equilibrium_dictionaries.append(self.model.mole_fraction_dict())

        return equilibrium_dictionaries

    def get_t50(self,
                species: str,
                criteria: Optional[str] = 'mass_frac',
                ):
        """
        Finds the half-life in seconds of the given species on either a mole fraction or mass fraction basis. Uses the
        initial conditions and reactor type when the class was initialized.

        Args:
            species (str): Cantera species name
            criteria (str): 'mol_frac' (50% reduction in X) or 'mass_frac' (50% reduction in Y)

        Returns:
            t50_list (list): List whose entries are the time [sec] of 50% conversion from initial amount.
        """
        t50_list = list()
        spc_index = self.spc_identifier_lookup[species]

        for condition in self.conditions:
            T0 = condition.T0.value_si
            try:
                V0 = self.conditions[0].V0.value_si
                P0 = None
            except AttributeError as e:
                P0 = condition.P0.value_si
                V0 = None
            self.reinitialize_simulation(T0=T0,
                                         P0=P0,
                                         X0=condition.mol_frac,
                                         V0=V0,
                                         )
            if criteria == 'mol_frac':
                x0 = self.model.X[spc_index]
            elif criteria == 'mass_frac':
                x0 = self.model.mass_fraction_dict()[species]
            else:
                raise ValueError(f'Invalid criteria: {criteria}')

            while True:
                self.cantera_simulation.step()

                if criteria == 'mol_frac':
                    x1 = self.model.X[spc_index]
                else:
                    x1 = self.model.mass_fraction_dict()[species]

                if x1 < x0 * 0.5:
                    break
            t50_list.append(self.cantera_simulation.time)

        return t50_list


register_simulate_adapter("CanteraConstantHP", CanteraConstantHP)
