"""
utility module for hashing gridfire reaction functions
"""
from __future__ import annotations
import typing
__all__: list[str] = ['mix_species', 'multiset_combine', 'splitmix64']
def mix_species(a: typing.SupportsInt, z: typing.SupportsInt) -> int:
    ...
def multiset_combine(acc: typing.SupportsInt, x: typing.SupportsInt) -> int:
    ...
def splitmix64(x: typing.SupportsInt) -> int:
    ...
