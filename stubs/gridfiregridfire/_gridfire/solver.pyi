"""
GridFire numerical solver bindings
"""
from __future__ import annotations
import collections.abc
import fourdst._phys.atomic
import gridfire._gridfire.engine
import gridfire._gridfire.engine.scratchpads
import gridfire._gridfire.type
import types
import typing
__all__: list[str] = ['GridSolver', 'GridSolverContext', 'MultiZoneDynamicNetworkSolver', 'PointSolver', 'PointSolverContext', 'PointSolverTimestepContext', 'SingleZoneDynamicNetworkSolver', 'SolverContextBase']
class GridSolver(MultiZoneDynamicNetworkSolver):
    def __init__(self, engine: gridfire._gridfire.engine.DynamicEngine, solver: SingleZoneDynamicNetworkSolver) -> None:
        """
        Initialize the GridSolver object.
        """
    def evaluate(self, solver_ctx: SolverContextBase, netIns: collections.abc.Sequence[gridfire._gridfire.type.NetIn]) -> list[gridfire._gridfire.type.NetOut]:
        """
        evaluate the dynamic engine using the dynamic engine class
        """
class GridSolverContext(SolverContextBase):
    detailed_logging: bool
    stdout_logging: bool
    zone_completion_logging: bool
    def __init__(self, ctx_template: gridfire._gridfire.engine.scratchpads.StateBlob) -> None:
        ...
    @typing.overload
    def clear_callback(self) -> None:
        ...
    @typing.overload
    def clear_callback(self, zone_idx: typing.SupportsInt) -> None:
        ...
    def init(self) -> None:
        ...
    def reset(self) -> None:
        ...
    @typing.overload
    def set_callback(self, callback: collections.abc.Callable[[...], None]) -> None:
        ...
    @typing.overload
    def set_callback(self, callback: collections.abc.Callable[[...], None], zone_idx: typing.SupportsInt) -> None:
        ...
class MultiZoneDynamicNetworkSolver:
    def evaluate(self, solver_ctx: SolverContextBase, netIns: collections.abc.Sequence[gridfire._gridfire.type.NetIn]) -> list[gridfire._gridfire.type.NetOut]:
        """
        evaluate the dynamic engine using the dynamic engine class for multiple zones (using openmp if available)
        """
class PointSolver(SingleZoneDynamicNetworkSolver):
    def __init__(self, engine: gridfire._gridfire.engine.DynamicEngine) -> None:
        """
        Initialize the PointSolver object.
        """
    def evaluate(self, solver_ctx: SolverContextBase, netIn: gridfire._gridfire.type.NetIn, display_trigger: bool = False, force_reinitialization: bool = False) -> gridfire._gridfire.type.NetOut:
        """
        evaluate the dynamic engine using the dynamic engine class
        """
class PointSolverContext:
    callback: collections.abc.Callable[[PointSolverTimestepContext], None] | None
    detailed_logging: bool
    stdout_logging: bool
    def __init__(self, engine_ctx: gridfire._gridfire.engine.scratchpads.StateBlob) -> None:
        ...
    def clear_context(self) -> None:
        ...
    def has_context(self) -> bool:
        ...
    def init(self) -> None:
        ...
    def init_context(self) -> None:
        ...
    def reset_all(self) -> None:
        ...
    def reset_cvode(self) -> None:
        ...
    def reset_user(self) -> None:
        ...
    @property
    def J(self) -> _generic_SUNMatrix:
        ...
    @property
    def LS(self) -> _generic_SUNLinearSolver:
        ...
    @property
    def Y(self) -> _generic_N_Vector:
        ...
    @property
    def YErr(self) -> _generic_N_Vector:
        ...
    @property
    def abs_tol(self) -> float:
        ...
    @abs_tol.setter
    def abs_tol(self, arg1: typing.SupportsFloat) -> None:
        ...
    @property
    def cvode_mem(self) -> types.CapsuleType:
        ...
    @property
    def engine_ctx(self) -> gridfire._gridfire.engine.scratchpads.StateBlob:
        ...
    @property
    def num_steps(self) -> int:
        ...
    @property
    def rel_tol(self) -> float:
        ...
    @rel_tol.setter
    def rel_tol(self, arg1: typing.SupportsFloat) -> None:
        ...
    @property
    def sun_ctx(self) -> SUNContext_:
        ...
class PointSolverTimestepContext:
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
    def state_ctx(self) -> gridfire._gridfire.engine.scratchpads.StateBlob:
        ...
    @property
    def t(self) -> float:
        ...
class SingleZoneDynamicNetworkSolver:
    def evaluate(self, solver_ctx: SolverContextBase, netIn: gridfire._gridfire.type.NetIn) -> gridfire._gridfire.type.NetOut:
        """
        evaluate the dynamic engine using the dynamic engine class for a single zone
        """
class SolverContextBase:
    pass
