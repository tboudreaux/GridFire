"""
GridFire plasma screening bindings
"""
from __future__ import annotations
import collections.abc
import fourdst._phys.atomic
import gridfire._gridfire.reaction
import typing
__all__: list[str] = ['BARE', 'BareScreeningModel', 'ScreeningModel', 'ScreeningType', 'WEAK', 'WeakScreeningModel', 'selectScreeningModel']
class BareScreeningModel:
    def __init__(self) -> None:
        ...
    def calculateScreeningFactors(self, reactions: gridfire._gridfire.reaction.ReactionSet, species: collections.abc.Sequence[fourdst._phys.atomic.Species], Y: collections.abc.Sequence[typing.SupportsFloat], T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> list[float]:
        """
        Calculate the bare plasma screening factors. This always returns 1.0 (bare)
        """
class ScreeningModel:
    pass
class ScreeningType:
    """
    Members:
    
      BARE
    
      WEAK
    """
    BARE: typing.ClassVar[ScreeningType]  # value = <ScreeningType.BARE: 0>
    WEAK: typing.ClassVar[ScreeningType]  # value = <ScreeningType.WEAK: 1>
    __members__: typing.ClassVar[dict[str, ScreeningType]]  # value = {'BARE': <ScreeningType.BARE: 0>, 'WEAK': <ScreeningType.WEAK: 1>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class WeakScreeningModel:
    def __init__(self) -> None:
        ...
    def calculateScreeningFactors(self, reactions: gridfire._gridfire.reaction.ReactionSet, species: collections.abc.Sequence[fourdst._phys.atomic.Species], Y: collections.abc.Sequence[typing.SupportsFloat], T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> list[float]:
        """
        Calculate the weak plasma screening factors using the Salpeter (1954) model.
        """
def selectScreeningModel(type: ScreeningType) -> ScreeningModel:
    """
    Select a screening model based on the specified type. Returns a pointer to the selected model.
    """
BARE: ScreeningType  # value = <ScreeningType.BARE: 0>
WEAK: ScreeningType  # value = <ScreeningType.WEAK: 1>
