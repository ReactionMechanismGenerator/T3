"""
RMS Simulator Adapter module
Used to run mechanism analysis with RMS at constant V
"""

import os
from typing import List, Optional, Type

from diffeqpy import de
from pyrms import rms

from t3.common import convert_termination_time_to_seconds
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter


class RMSConstantP(SimulateAdapter):
    """
    RMSConstantP is an adapter for the abstract class SimulateAdapter that runs homogeneous gas reactions
    under isobaric conditions.

    Args:
        t3 (dict): The T3.t3 attribute, which is a dictionary containing the t3 block from the input yaml or API.
        rmg (dict): The T3.rmg attribute, which is a dictionary containing the rmg block from the input yaml or API.
        paths (dict): The T3.paths attribute, which is a dictionary containing relevant paths.
        logger (Logger): Instance of T3's Logger class.
        atol (float): The absolute tolerance used when integrating.
        rtol (float): The relative tolerance used when integrating.
        observable_list (Optional[list]): Species used for SA. Entries are species labels as strings. Example: ['OH']
        sa_atol (Optional[float]): The absolute tolerance used when performing sensitivity analysis.
        sa_atol (Optional[float]): The relative tolerance used when performing sensitivity analysis.
        global_observables (Optional[List[str]]): List of global observables ['IgD', 'ESR', 'SL'] used by Cantera adapters.

    Attributes:
        atol (float): The absolute tolerance used when integrating.
        bsol (Simulation object): Simulation object from the RMS Simulation function.
        global_observables (Optional[List[str]]): List of global observables ['IgD', 'ESR', 'SL'] used by Cantera adapters.
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
                 sa_atol: Optional[float] = 1e-6,
                 sa_rtol: Optional[float] = 1e-4,
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
        self.rms_directory = os.path.join(self.paths['RMG'], 'rms')
        self.rms_file = None
        self.phaseDict = None
        self.spcs = None
        self.rxns = None
        self.initialconds = dict()
        self.phase = None
        self.sol = None
        self.bsol = None

        self.set_up()


    def set_up(self):
        """
        Read in the rms yaml file, obtain simulation conditions from the T3 input file, and simulate the job.
        """
        if len(self.observable_list):
            self.logger.info('Running a simulation with SA using RMS...')
        else:
            self.logger.info('Running a simulation using RMS...')

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
        domain, y0, p = rms.ConstantPDomain(phase=self.phase, initialconds=self.initialconds)

        solver = de.CVODE_BDF()
        t_final = convert_termination_time_to_seconds(self.rmg['reactors'][0]['termination_time'])

        if len(self.observable_list):
            react = rms.Reactor(domain, y0, (0.0, t_final), forwardsensitivities=True, p=p)
            atol = [self.atol]*domain.indexes[-1]  # set atol for the variables being simulated
            total_variables = domain.indexes[-1] + len(p)*domain.indexes[-1]
            atol.extend([self.sa_atol]*(total_variables - domain.indexes[-1]))  # set atol for SA
            self.sol = de.solve(react.ode, solver, abstol=atol, reltol=self.rtol)
        else:
            react = rms.Reactor(domain, y0, (0.0, t_final), forwardsensitivities=False, p=p)
            self.sol = de.solve(react.ode, solver, abstol=self.atol, reltol=self.rtol)
        self.bsol = rms.Simulation(self.sol, domain)


    def get_sa_coefficients(self):
        """
        Obtain the SA coefficients.

        Returns:
             sa_dict: a SA dictionary, whose structure is given in the docstring for main.py
        """
        pass
        # sa_dict = {'thermo': dict(), 'kinetics': dict(), 'time': list(self.sol.t)}
        # species_labels = [spc.name for spc in self.spcs]
        # num_rxns = len(self.rxns)
        # for observable in self.observable_list:
        #     # get thermo SA coefficients
        #     sa_dict['thermo'][observable] = dict()
        #     for label in species_labels:
        #         sa_dict['thermo'][observable][label] = list()
        #         for t in sa_dict['time']:
        #             sa_dict['thermo'][observable][label].append(rms.getconcentrationsensitivity(
        #                 self.bsol, observable, label, t))
        #
        #     # get kinetics SA coefficients
        #     sa_dict['kinetics'][observable] = dict()
        #     for rxn_number in range(num_rxns):
        #         sa_dict['kinetics'][observable][label] = list()
        #         for t in sa_dict['time']:
        #             sa_dict['kinetics'][observable][label].append(rms.getconcentrationsensitivity(
        #                 self.bsol, observable, rxn_number, t))
        # return sa_dict


register_simulate_adapter("RMSConstantP", RMSConstantP)
