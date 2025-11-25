"""
GridFire numerical solver bindings
"""
from __future__ import annotations
import collections.abc
import fourdst._phys.atomic
import gridfire._gridfire.engine
import gridfire._gridfire.type
import typing
__all__: list[str] = ['CVODESolverStrategy', 'CVODETimestepContext', 'DynamicNetworkSolverStrategy', 'SolverContextBase']
class CVODESolverStrategy(DynamicNetworkSolverStrategy):
    def __init__(self, engine: gridfire._gridfire.engine.DynamicEngine) -> None:
        """
        Initialize the CVODESolverStrategy object.
        """
    def evaluate(self, netIn: gridfire._gridfire.type.NetIn, display_trigger: bool = False) -> gridfire._gridfire.type.NetOut:
        """
        evaluate the dynamic engine using the dynamic engine class
        """
    def get_absTol(self) -> float:
        """
        Get the absolute tolerance for the CVODE solver.
        """
    def get_relTol(self) -> float:
        """
        Get the relative tolerance for the CVODE solver.
        """
    def get_stdout_logging_enabled(self) -> bool:
        """
        Check if solver logging to standard output is enabled.
        """
    def set_absTol(self, absTol: typing.SupportsFloat) -> None:
        """
        Set the absolute tolerance for the CVODE solver.
        """
    def set_callback(self, cb: collections.abc.Callable[[CVODETimestepContext], None]) -> None:
        """
        Set a callback function which will run at the end of every successful timestep
        """
    def set_relTol(self, relTol: typing.SupportsFloat) -> None:
        """
        Set the relative tolerance for the CVODE solver.
        """
    def set_stdout_logging_enabled(self, logging_enabled: bool) -> None:
        """
        Enable logging to standard output.
        """
class CVODETimestepContext(SolverContextBase):
    @property
    def T9(self) -> float:
        ...
    @property
    def currentConvergenceFailures(self) -> int:
        ...
    @property
    def currentNonlinearIterations(self) -> int:
        ...
    @property
    def dt(self) -> float:
        ...
    @property
    def engine(self) -> gridfire._gridfire.engine.DynamicEngine:
        ...
    @property
    def last_step_time(self) -> float:
        ...
    @property
    def networkSpecies(self) -> list[fourdst._phys.atomic.Species]:
        ...
    @property
    def num_steps(self) -> int:
        ...
    @property
    def rho(self) -> float:
        ...
    @property
    def state(self) -> list[float]:
        ...
    @property
    def t(self) -> float:
        ...
class DynamicNetworkSolverStrategy:
    def describe_callback_context(self) -> list[tuple[str, str]]:
        """
        Get a structure representing what data is in the callback context in a human readable format
        """
    def evaluate(self, netIn: gridfire._gridfire.type.NetIn) -> gridfire._gridfire.type.NetOut:
        """
        evaluate the dynamic engine using the dynamic engine class
        """
class SolverContextBase:
    pass
