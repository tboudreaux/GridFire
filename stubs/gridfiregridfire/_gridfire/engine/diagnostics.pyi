"""
A submodule for engine diagnostics
"""
from __future__ import annotations
import collections.abc
import fourdst._phys.composition
import gridfire._gridfire.engine
import gridfire._gridfire.engine.scratchpads
import typing
__all__: list[str] = ['inspect_jacobian_stiffness', 'inspect_species_balance', 'report_limiting_species']
def inspect_jacobian_stiffness(ctx: gridfire._gridfire.engine.scratchpads.StateBlob, engine: gridfire._gridfire.engine.DynamicEngine, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, json: bool) -> ... | None:
    ...
def inspect_species_balance(ctx: gridfire._gridfire.engine.scratchpads.StateBlob, engine: gridfire._gridfire.engine.DynamicEngine, species_name: str, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, json: bool) -> ... | None:
    ...
def report_limiting_species(ctx: gridfire._gridfire.engine.scratchpads.StateBlob, engine: gridfire._gridfire.engine.DynamicEngine, Y_full: collections.abc.Sequence[typing.SupportsFloat], E_full: collections.abc.Sequence[typing.SupportsFloat], relTol: typing.SupportsFloat, absTol: typing.SupportsFloat, top_n: typing.SupportsInt, json: bool) -> ... | None:
    ...
