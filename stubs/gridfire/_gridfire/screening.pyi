"""
GridFire plasma screening bindings
"""
from __future__ import annotations
import collections.abc
import fourdst._phys.atomic
import gridfire._gridfire.reaction
import typing
__all__: list[str] = ['BARE', 'BLENDED', 'BareScreeningModel', 'BlendedScreeningModel', 'CD09', 'CD09ScreeningModel', 'ScreeningBlendWindow', 'ScreeningModel', 'ScreeningType', 'WEAK', 'WeakScreeningModel', 'selectScreeningModel']
class BareScreeningModel:
    def __init__(self) -> None:
        ...
    def calculateScreeningFactors(self, reactions: gridfire._gridfire.reaction.ReactionSet, species: collections.abc.Sequence[fourdst._phys.atomic.Species], Y: collections.abc.Sequence[typing.SupportsFloat], T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> list[float]:
        """
        Calculate the bare plasma screening factors. This always returns 1.0 (bare)
        """
class BlendedScreeningModel:
    def __init__(self, window: ScreeningBlendWindow = ...) -> None:
        """
        Salpeter weak screening handing over to Chugunov & DeWitt (2009) across the blend window
        """
    def calculateScreeningFactors(self, reactions: gridfire._gridfire.reaction.ReactionSet, species: collections.abc.Sequence[fourdst._phys.atomic.Species], Y: collections.abc.Sequence[typing.SupportsFloat], T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> list[float]:
        """
        Calculate plasma screening factors, blending Salpeter and Chugunov & DeWitt (2009).
        """
    def weight(self, H: typing.SupportsFloat) -> float:
        """
        Blend weight for a weak-screening exponent H (0: weak only, 1: CD09 only)
        """
class CD09ScreeningModel:
    def __init__(self) -> None:
        ...
    def calculateScreeningFactors(self, reactions: gridfire._gridfire.reaction.ReactionSet, species: collections.abc.Sequence[fourdst._phys.atomic.Species], Y: collections.abc.Sequence[typing.SupportsFloat], T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> list[float]:
        """
        Calculate plasma screening factors with Chugunov & DeWitt (2009), valid from weak to strong coupling.
        """
class ScreeningBlendWindow:
    strongExponent: float
    weakExponent: float
    def __init__(self) -> None:
        ...
class ScreeningModel:
    pass
class ScreeningType:
    """
    Members:
    
      BARE
    
      WEAK
    
      CD09
    
      BLENDED
    """
    BARE: typing.ClassVar[ScreeningType]  # value = <ScreeningType.BARE: 0>
    WEAK: typing.ClassVar[ScreeningType]  # value = <ScreeningType.WEAK: 1>
    CD09: typing.ClassVar[ScreeningType]  # value = <ScreeningType.CD09: 2>
    BLENDED: typing.ClassVar[ScreeningType]  # value = <ScreeningType.BLENDED: 3>
    __members__: typing.ClassVar[dict[str, ScreeningType]]  # value = {'BARE': <ScreeningType.BARE: 0>, 'WEAK': <ScreeningType.WEAK: 1>, 'CD09': <ScreeningType.CD09: 2>, 'BLENDED': <ScreeningType.BLENDED: 3>}
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
def selectScreeningModel(type: ScreeningType, window: ScreeningBlendWindow = ...) -> ScreeningModel:
    """
    Select a screening model based on the specified type. Returns a pointer to the selected model.
    """
BARE: ScreeningType  # value = <ScreeningType.BARE: 0>
WEAK: ScreeningType  # value = <ScreeningType.WEAK: 1>
CD09: ScreeningType  # value = <ScreeningType.CD09: 2>
BLENDED: ScreeningType  # value = <ScreeningType.BLENDED: 3>
