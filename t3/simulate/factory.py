"""
A module for generating simulate adapters.
"""

from typing import TYPE_CHECKING, List, Optional, Type

from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
if TYPE_CHECKING:
    from t3.simulate.adapter import SimulateAdapter

_registered_simulate_adapters = {}


def register_simulate_adapter(simulator: str,
                              simulate_class: SimulateAdapter,
                              ) -> None:
    """
    A register for the simulate adapters.

    Args:
        simulator (str): A string representation for a simulate adapter.
        simulate_class (Type[SimulateAdapter]): The simulate adapter class.

    Raises:
        TypeError: If ``simulate_class`` is not a ``SimulateAdapter`` instance.
    """
    if not issubclass(simulate_class, SimulateAdapter):
        raise TypeError(f'simulate_class is not a SimulateAdapter, got {simulate_class} which is a {type(simulate_class)}')
    _registered_simulate_adapters[simulator] = simulate_class


def simulate_factory(simulate_method: str,
                     t3: dict,
                     rmg: dict,
                     paths: dict,
                     logger: Type[Logger],
                     atol: float,
                     rtol: float,
                     observable_list: Optional[list] = None,
                     sa_atol: float = 1e-6,
                     sa_rtol: float = 1e-4,
                     global_observables: Optional[List[str]] = None,
                     ) -> SimulateAdapter:
    """
    A factory generating the simulate adapter corresponding to ``simulate_adapter``.

    Args:
        simulate_method (str): The simulate adapter name. Examples: 'RMGConstantTP', 'RMSConstantTP',
                               or 'CanteraConstantTP'.
        t3 (dict): The T3.t3 attribute, which is a dictionary containing the t3 block from the input yaml or API.
        rmg (dict): The T3.rmg attribute, which is a dictionary containing the rmg block from the input yaml or API.
        paths (dict): The T3.paths attribute, which is a dictionary containing relevant paths.
        logger (Type[Logger]): The current T3 Logger instance.
        atol (float): The absolute tolerance used when integrating.
        rtol (float): The relative tolerance used when integrating.
        observable_list (Optional[list]): Species used for SA. Entries are species labels as strings. Example: ['OH']
        sa_atol (float, optional): The absolute tolerance used when performing sensitivity analysis.
        sa_rtol (float, optional): The relative tolerance used when performing sensitivity analysis.
        global_observables (Optional[List[str]]): List of global observables ['IgD', 'ESR', 'SL'] used by Cantera.

    Raises:
        ValueError: If the provided simulate_method is not in the keys for the _registered_simulate_adapters dictionary.

    Returns:
        Type[SimulateAdapter]: The requested SimulateAdapter child, initialized with the respective arguments.
    """

    if simulate_method not in _registered_simulate_adapters.keys():
        raise ValueError(f'The "simulate_method" argument of {simulate_method} was not present in the keys for the '
                         f'_registered_simulate_adapters dictionary: {list(_registered_simulate_adapters.keys())}'
                         f'\nPlease check that the simulate adapter was registered properly.')

    simulate_adapter_class = _registered_simulate_adapters[simulate_method](t3=t3,
                                                                            rmg=rmg,
                                                                            paths=paths,
                                                                            logger=logger,
                                                                            atol=atol,
                                                                            rtol=rtol,
                                                                            observable_list=observable_list,
                                                                            sa_atol=sa_atol,
                                                                            sa_rtol=sa_rtol,
                                                                            global_observables=global_observables,
                                                                            )
    return simulate_adapter_class
