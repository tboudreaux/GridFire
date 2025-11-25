"""
GridFire utility method bindings
"""
from __future__ import annotations
import fourdst._phys.composition
import gridfire._gridfire.engine
import typing
from . import hashing
__all__: list[str] = ['formatNuclearTimescaleLogString', 'hash_atomic', 'hash_reaction', 'hashing']
def formatNuclearTimescaleLogString(engine: gridfire._gridfire.engine.DynamicEngine, Y: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> str:
    """
    Format a string for logging nuclear timescales based on temperature, density, and energy generation rate.
    """
def hash_atomic(a: typing.SupportsInt, z: typing.SupportsInt) -> int:
    ...
def hash_reaction(reaction: ...) -> int:
    ...
