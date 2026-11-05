"""
GridFire configuration bindings
"""
from __future__ import annotations
import gridfire._gridfire.engine
import gridfire._gridfire.screening
import typing
__all__: list[str] = ['AdaptiveEngineViewConfig', 'BoundaryFluxConfig', 'EngineConfig', 'EngineViewConfig', 'GridFireConfig', 'PointSolverConfig', 'SolverConfig', 'TriggerConfig']
class AdaptiveEngineViewConfig:
    def __init__(self) -> None:
        ...
    @property
    def relativeCullingThreshold(self) -> float:
        ...
    @relativeCullingThreshold.setter
    def relativeCullingThreshold(self, arg0: typing.SupportsFloat) -> None:
        ...
class BoundaryFluxConfig:
    def __init__(self) -> None:
        ...
    @property
    def absoluteThreshold(self) -> float:
        ...
    @absoluteThreshold.setter
    def absoluteThreshold(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def relativeThreshold(self) -> float:
        ...
    @relativeThreshold.setter
    def relativeThreshold(self, arg0: typing.SupportsFloat) -> None:
        ...
class EngineConfig:
    jacobian: gridfire._gridfire.engine.JacobianMethod
    screening: gridfire._gridfire.screening.ScreeningType
    screeningBlend: gridfire._gridfire.screening.ScreeningBlendWindow
    passiveSpecies: list[str]
    views: EngineViewConfig
    def __init__(self) -> None:
        ...
class EngineViewConfig:
    adaptiveEngineView: AdaptiveEngineViewConfig
    def __init__(self) -> None:
        ...
class GridFireConfig:
    engine: EngineConfig
    solver: SolverConfig
    def __init__(self) -> None:
        ...
class PointSolverConfig:
    clampNegativeIterates: bool
    trigger: TriggerConfig
    def __init__(self) -> None:
        ...
    @property
    def absTol(self) -> float:
        ...
    @absTol.setter
    def absTol(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def maxStepFailureRestarts(self) -> int:
        ...
    @maxStepFailureRestarts.setter
    def maxStepFailureRestarts(self, arg0: typing.SupportsInt) -> None:
        ...
    @property
    def relTol(self) -> float:
        ...
    @relTol.setter
    def relTol(self, arg0: typing.SupportsFloat) -> None:
        ...
class SolverConfig:
    pointSolver: PointSolverConfig
    def __init__(self) -> None:
        ...
class TriggerConfig:
    boundaryFlux: BoundaryFluxConfig
    def __init__(self) -> None:
        ...
    @property
    def maxConvergenceFailures(self) -> float:
        ...
    @maxConvergenceFailures.setter
    def maxConvergenceFailures(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def offDiagonalThreshold(self) -> float:
        ...
    @offDiagonalThreshold.setter
    def offDiagonalThreshold(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def timestepCollapseRatio(self) -> float:
        ...
    @timestepCollapseRatio.setter
    def timestepCollapseRatio(self, arg0: typing.SupportsFloat) -> None:
        ...
