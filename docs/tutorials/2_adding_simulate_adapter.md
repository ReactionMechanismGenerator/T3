# 2. Adding a Simulate Adapter

This tutorial demonstrates how to create a new simulate adapter for T3 to use when simulating
a mechanism. The currently implemented Cantera, RMG, and RMS adapters provide good examples for this process, 
though this tutorial explains the steps in more detail by walking through each piece of the code block below. 
First, create a file under `T3/t3/simulate/` where we will write the new simulate adapter class.
Next, add the required imports as shown below. Additional packages may be needed depending on the task, but at
a minimum, these are required.

```Python hl_lines="1 2 3 4"
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter

class NewSimulator(SimulateAdapter):
    """
    Briefly summarize what conditions and/or assumptions are used in this simulator.

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
        <varies by class> 
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
        
        # initialize attributes
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
        
        # run any additional setup required by the class
        self.set_up()


    def set_up(self):
        """
        Set up the class, possibly by reading an input file, initializing values, and/or setting up the domain.
        """
        <Implement code here>    


    def simulate(self):
        """
        Simulate a job to obtain species profiles. Run sensitivity analysis if requested.
        """
        <Implement code here>


    def get_sa_coefficients(self):
        """
        Obtain the sensitivity analysis coefficients.
        
        Returns:
             sa_dict (dict): a SA dictionary, whose structure is given in the docstring for T3/t3/main.py
        """
        <Implement code here>


    def get_idt_by_T(self):
        """
        Finds the ignition point by approximating dT/dt as a first order forward difference
        and then finds the point of maximum slope.

        Returns:
            idt_dict (dict): Dictionary whose keys include 'idt' and 'idt_index' and whose values are lists of
                             the ignition delay time in seconds and index at which idt occurs respectively.
        """
        <Implement code here>



register_simulate_adapter("NewSimulator", NewSimulator)

```

The new class must inherit from the abstract adapter class in `T3/t3/simulate/adapter.py`. All simulate adapters
accept the same arguments, allowing them to be used in a standardized way. Thus, the `Args` section in the docstring
and the `__init__` method will be constant for all adapters.


```Python hl_lines="5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53"
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter

class NewSimulator(SimulateAdapter):
    """
    Briefly summarize what conditions and/or assumptions are used in this simulator.

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
        <varies by class> 
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
        
        # initialize attributes
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
        
        # run any additional setup required by the class
        self.set_up()


    def set_up(self):
        """
        Set up the class, possibly by reading an input file, initializing values, and/or setting up the domain.
        """
        <Implement code here>    


    def simulate(self):
        """
        Simulate a job to obtain species profiles. Run sensitivity analysis if requested.
        """
        <Implement code here>


    def get_sa_coefficients(self):
        """
        Obtain the sensitivity analysis coefficients.
        
        Returns:
             sa_dict (dict): a SA dictionary, whose structure is given in the docstring for T3/t3/main.py
        """
        <Implement code here>


    def get_idt_by_T(self):
        """
        Finds the ignition point by approximating dT/dt as a first order forward difference
        and then finds the point of maximum slope.

        Returns:
            idt_dict (dict): Dictionary whose keys include 'idt' and 'idt_index' and whose values are lists of
                             the ignition delay time in seconds and index at which idt occurs respectively.
        """
        <Implement code here>


register_simulate_adapter("NewSimulator", NewSimulator)

```


The next step is to implement the following methods that were inherited by the abstract class: `set_up()`,
`simulate()`, `get_sa_coefficients()`, and  `get_idt_by_T()`. These methods are present in all simulators, and 
 if they return anything, do so in a standard format. This standardizes the use of different simulators, making it 
 easy for the user to simulate their mechanism with different adapters. 

```Python hl_lines="54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86"
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter

class NewSimulator(SimulateAdapter):
    """
    Briefly summarize what conditions and/or assumptions are used in this simulator.

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
        <varies by class> 
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
        
        # initialize attributes
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
        
        # run any additional setup required by the class
        self.set_up()


    def set_up(self):
        """
        Set up the class, possibly by reading an input file, initializing values, and/or setting up the domain.
        """
        <Implement code here>    


    def simulate(self):
        """
        Simulate a job to obtain species profiles. Run sensitivity analysis if requested.
        """
        <Implement code here>


    def get_sa_coefficients(self):
        """
        Obtain the sensitivity analysis coefficients.
        
        Returns:
             sa_dict (dict): a SA dictionary, whose structure is given in the docstring for T3/t3/main.py
        """
        <Implement code here>


    def get_idt_by_T(self):
        """
        Finds the ignition point by approximating dT/dt as a first order forward difference
        and then finds the point of maximum slope.

        Returns:
            idt_dict (dict): Dictionary whose keys include 'idt' and 'idt_index' and whose values are lists of
                             the ignition delay time in seconds and index at which idt occurs respectively.
        """
        <Implement code here>



register_simulate_adapter("NewSimulator", NewSimulator)

```

Finally, register the adapter with the factory at the bottom of the file. 
Also, initialize the simulator by importing it in `T3/t3/simulate/__init__.py`.

```Python hl_lines="87 88 89"
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter

class NewSimulator(SimulateAdapter):
    """
    Briefly summarize what conditions and/or assumptions are used in this simulator.

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
        <varies by class> 
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
        
        # initialize attributes
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
        
        # run any additional setup required by the class
        self.set_up()


    def set_up(self):
        """
        Set up the class, possibly by reading an input file, initializing values, and/or setting up the domain.
        """
        <Implement code here>    


    def simulate(self):
        """
        Simulate a job to obtain species profiles. Run sensitivity analysis if requested.
        """
        <Implement code here>


    def get_sa_coefficients(self):
        """
        Obtain the sensitivity analysis coefficients.
        
        Returns:
             sa_dict (dict): a SA dictionary, whose structure is given in the docstring for T3/t3/main.py
        """
        <Implement code here>


    def get_idt_by_T(self):
        """
        Finds the ignition point by approximating dT/dt as a first order forward difference
        and then finds the point of maximum slope.

        Returns:
            idt_dict (dict): Dictionary whose keys include 'idt' and 'idt_index' and whose values are lists of
                             the ignition delay time in seconds and index at which idt occurs respectively.
        """
        <Implement code here>



register_simulate_adapter("NewSimulator", NewSimulator)

```
