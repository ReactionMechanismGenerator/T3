"""
RMS Simulator Adapter module
Used to run mechanism analysis with RMS at constant TP
"""

import os
from typing import List, Optional, Type

from diffeqpy import de
from pyrms import rms

from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter


class RMSConstantTP(SimulateAdapter):
    """
    RMSConstantTP is an adapter for the abstract class SimulateAdapter that runs homogeneous gas reactions
    under isothermal and isobaric conditions.

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
        atol (float): The absolute tolerance used when integrating during an RMG iteration.
        bsol (Simulation object): Simulation object from the RMS Simulation function.
        global_observables (Optional[List[str]]): List of global observables ['IgD', 'ESR', 'SL'] used by Cantera adapters.
        initialconds (dict): Dictionary containing initial conditions for the simulation.
        logger (Logger): Instance of T3's Logger class.
        observable_list (list): Species used for SA. Entries are species labels as strings. Example: ['OH']
        paths (dict): The T3.paths attribute, which is a dictionary containing relevant paths.
        phase (object): The RMS phase, such as IdealGas or IdealDiluteSolution.
        phaseDict (dict): Result of reading in the rms yaml file.
        rmg (dict): The T3.rmg attribute, which is a dictionary containing the rmg block from the input yaml or API.
        rms_directory (str): Path to the "RMG/rms" directory for the current iteration.
        rms_file (str): Path to the rms yaml for the current iteration.
        rtol (float): The relative tolerance used when integrating during an RMG iteration.
        rxns (list): List of AbstractReaction object.
        sa_atol (float): The absolute tolerance used when performing sensitivity analysis.
        sa_rtol (float): The relative tolerance used when performing sensitivity analysis.
        sol (ODESolution object): Output from the DifferentialEquations package.
        spcs (list): List of AbstractSpecies object.
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
        self.phase = rms.IdealGas(self.spcs, self.rxns, name="phase")

        # Define the domain (encodes how system thermodynamic properties calculated)
        domain, y0 = rms.ConstantTPDomain(phase=self.phase, initialconds=self.initialconds,
                                          sensitivity=True if len(self.observable_list) else False)

        solver = de.CVODE_BDF()
        react = rms.Reactor(domain, y0, (0.0, self.rmg['reactors'][0]['termination_time']))

        self.sol = de.solve(react.ode, solver, abstol=self.atol, reltol=self.rtol)
        self.bsol = rms.Simulation(self.sol, domain)


    def get_sa_coefficients(self):
        """
        Obtain the SA coefficients.

        Returns:
             sa_dict: a SA dictionary, whose structure is given in the docstring for main.py
        """
        pass


register_simulate_adapter("RMSConstantTP", RMSConstantTP)
