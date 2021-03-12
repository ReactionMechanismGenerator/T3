"""
RMS Simulator Adapter module
Used to run mechanism analysis with RMS that simulates ideal gases at constant T and constant P
"""

import os
from typing import List, Optional, Type

from t3.common import convert_termination_time_to_seconds
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter

HAS_RMS = True
try:
    from diffeqpy import de
    from pyrms import rms
except (ImportError, ModuleNotFoundError):
    # diffeqpy and/or pyrms are missing
    HAS_RMS = False


class RMSConstantTP(SimulateAdapter):
    """
    RMSConstantTP is an adapter for the abstract class SimulateAdapter that simulates ideal gases
    in a batch reactor at constant temperature and constant pressure.

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
        atol (float): The absolute tolerance used when integrating.
        bsol (Simulation object): Simulation object from the RMS Simulation function.
        global_observables (List[str]): List of global observables ['IgD', 'ESR', 'SL'] used by Cantera adapters.
        initialconds (dict): Dictionary containing initial conditions for the simulation.
        logger (Logger): Instance of T3's Logger class.
        observable_list (list): Species used for SA. Entries are species labels as strings. Example: ['OH']
        paths (dict): The T3.paths attribute, which is a dictionary containing relevant paths.
        phase (object): The RMS phase, such as IdealGas or IdealDiluteSolution.
        phaseDict (dict): Result of reading in the rms yaml file, which is a dictionary with one key: 'phase'
        rmg (dict): The T3.rmg attribute, which is a dictionary containing the rmg block from the input yaml or API.
        rms_directory (str): Path to the "RMG/rms" directory for the current iteration.
        rms_file (str): Path to the rms yaml for the current iteration.
        rtol (float): The relative tolerance used when integrating.
        rxns (list): List of AbstractReaction object. Each entry has attributes, such asindex, reactants, reactantinds,
                     products, productsinds, radicalchange, pairs.
        sa_atol (float): The absolute tolerance used when performing sensitivity analysis.
        sa_rtol (float): The relative tolerance used when performing sensitivity analysis.
                         de.solve() only accepts scalar values for reltol so this value is not used when doing SA.
                         In contrast, abstol can accept a vector of values, which allows different tolerances
                         to be used for simulating and for performing SA.
        sol (ODESolution object): Output from the DifferentialEquations package.
        spcs (list): List of AbstractSpecies object. Each entry has attributes, such as name, index, inchi, smiles,
                     thermo, atomnums, bondnum, diffusion, radius, radicalelectrons.
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
        self.constantspecies = [species['label'] for species in self.rmg['species'] if species['constant']]

        # cantera argument not used by this adapter
        self.global_observables = global_observables

        # initialize other attributes
        self.rms_directory = self.paths['RMS']
        self.rms_file = None
        self.phaseDict = None
        self.spcs = None
        self.rxns = None
        self.initialconds = dict()
        self.phase = None
        self.domain = None
        self.y0 = None
        self.p = None
        self.sol = None
        self.bsol = None

        self.set_up()

    def set_up(self):
        """
        Read in the rms yaml file, obtain simulation conditions from the T3 input file, and set up the domain.
        """

        # read in rms yaml file to obtain species and reactions
        max_filenum = -1
        for file in os.listdir(self.rms_directory):
            if file.endswith(".rms"):
                filename = file.split('.')[0]
                filenum = int(filename.split('chem')[-1])
                if filenum > max_filenum:
                    max_filenum = filenum

        self.rms_file = os.path.join(self.rms_directory, "chem" + str(max_filenum) + ".rms")
        self.phaseDict = rms.readinput(self.rms_file)
        self.spcs = self.phaseDict["phase"]["Species"]
        self.rxns = self.phaseDict["phase"]["Reactions"]

        # obtain initial mole fraction values from the T3 input file or API
        for input_species in self.rmg['species']:
            self.initialconds.update({input_species['label']: input_species['concentration']})

        # obtain initial T, P
        self.initialconds.update({"P": self.rmg['reactors'][0]['P'], "T": self.rmg['reactors'][0]['T']})

        # Define the phase (how species thermodynamic and kinetic properties calculated)
        self.phase = rms.IdealGas(self.spcs, self.rxns, name="gas")

        # Define the domain (encodes how system thermodynamic properties calculated)
        self.domain, self.y0, self.p = rms.ConstantTPDomain(phase=self.phase,
                                                            initialconds=self.initialconds,
                                                            constantspecies=self.constantspecies,
                                                            )

    def simulate(self):
        """
        Simulate the mechanism with RMS
        """
        if not HAS_RMS:
            self.logger.error('Missing diffeqpy and/or pyrms. Please ensure these dependencies were'
                              'correctly installed.')

        if len(self.observable_list):
            self.logger.info('Running a simulation with SA using RMSConstantTP...')
        else:
            self.logger.info('Running a simulation using RMSConstantTP...')

        solver = de.CVODE_BDF()
        t_final = convert_termination_time_to_seconds(self.rmg['reactors'][0]['termination_time'])

        if len(self.observable_list):
            react = rms.Reactor(self.domain, self.y0, (0.0, t_final), forwardsensitivities=True, p=self.p)
            atol = [self.atol] * self.domain.indexes[-1]  # set atol for the variables being simulated
            total_variables = self.domain.indexes[-1] + len(self.p) * self.domain.indexes[-1]
            atol.extend([self.sa_atol] * (total_variables - self.domain.indexes[-1]))  # set atol for SA
            self.sol = de.solve(react.ode, solver, abstol=atol, reltol=self.rtol)
        else:
            react = rms.Reactor(self.domain, self.y0, (0.0, t_final), forwardsensitivities=False, p=self.p)
            self.sol = de.solve(react.ode, solver, abstol=self.atol, reltol=self.rtol)
        self.bsol = rms.Simulation(self.sol, self.domain)

    def get_sa_coefficients(self):
        """
        Obtain the SA coefficients.

        Returns:
             sa_dict (dict): a SA dictionary, whose structure is given in the docstring for T3/t3/main.py
        """
        sa_dict = {'thermo': dict(), 'kinetics': dict(), 'time': list(self.sol.t)}
        species_labels = [spc.name for spc in self.spcs]
        num_rxns = len(self.rxns)
        for observable in self.observable_list:
            # get thermo SA coefficients
            sa_dict['thermo'][observable] = dict()
            for label in species_labels:
                sa_dict['thermo'][observable][label] = list()
                for t in sa_dict['time']:
                    sa_dict['thermo'][observable][label].append(rms.getconcentrationsensitivity(
                        self.bsol, observable, label, t))

            # get kinetics SA coefficients
            sa_dict['kinetics'][observable] = dict()
            for rxn_number in range(num_rxns):
                sa_dict['kinetics'][observable][rxn_number] = list()
                for t in sa_dict['time']:
                    sa_dict['kinetics'][observable][rxn_number].append(rms.getconcentrationsensitivity(
                        self.bsol, observable, rxn_number, t))
        return sa_dict

    def get_idt_by_T(self):
        """
        Finds the ignition point by approximating dT/dt as a first order forward difference
        and then finds the point of maximum slope. Since this adapter simulates at constant T, this method
        returns a dictionary whose values are empty lists.

        Returns:
            idt_dict (dict): Dictionary whose keys include 'idt' and 'idt_index' and whose values are lists of
                             the ignition delay time in seconds and index at which idt occurs respectively.
        """
        idt_dict = {'idt': list(),
                    'idt_index': list(),
                    }

        return idt_dict


register_simulate_adapter("RMSConstantTP", RMSConstantTP)
