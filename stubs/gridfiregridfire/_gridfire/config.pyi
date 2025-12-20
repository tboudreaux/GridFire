"""
GridFire configuration bindings
"""
from __future__ import annotations
import typing
__all__: list[str] = ['AdaptiveEngineViewConfig', 'CVODESolverConfig', 'EngineConfig', 'EngineViewConfig', 'GridFireConfig', 'SolverConfig']
class AdaptiveEngineViewConfig:
    def __init__(self) -> None:
        ...
    @property
    def relativeCullingThreshold(self) -> float:
        ...
    @relativeCullingThreshold.setter
    def relativeCullingThreshold(self, arg0: typing.SupportsFloat) -> None:
        ...
class CVODESolverConfig:
    def __init__(self) -> None:
        ...
    @property
    def absTol(self) -> float:
        ...
    @absTol.setter
    def absTol(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def relTol(self) -> float:
        ...
    @relTol.setter
    def relTol(self, arg0: typing.SupportsFloat) -> None:
        ...
class EngineConfig:
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
class SolverConfig:
    cvode: CVODESolverConfig
    def __init__(self) -> None:
        ...
