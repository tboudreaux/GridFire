"""
GridFire partition function bindings
"""
from __future__ import annotations
import collections.abc
import typing
__all__: list[str] = ['BasePartitionType', 'CompositePartitionFunction', 'GroundState', 'GroundStatePartitionFunction', 'PartitionFunction', 'RauscherThielemann', 'RauscherThielemannPartitionDataRecord', 'RauscherThielemannPartitionFunction', 'basePartitionTypeToString', 'stringToBasePartitionType']
class BasePartitionType:
    """
    Members:
    
      RauscherThielemann
    
      GroundState
    """
    GroundState: typing.ClassVar[BasePartitionType]  # value = <BasePartitionType.GroundState: 1>
    RauscherThielemann: typing.ClassVar[BasePartitionType]  # value = <BasePartitionType.RauscherThielemann: 0>
    __members__: typing.ClassVar[dict[str, BasePartitionType]]  # value = {'RauscherThielemann': <BasePartitionType.RauscherThielemann: 0>, 'GroundState': <BasePartitionType.GroundState: 1>}
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
class CompositePartitionFunction:
    @typing.overload
    def __init__(self, partitionFunctions: collections.abc.Sequence[BasePartitionType]) -> None:
        """
        Create a composite partition function from a list of base partition types.
        """
    @typing.overload
    def __init__(self, arg0: CompositePartitionFunction) -> None:
        """
        Copy constructor for CompositePartitionFunction.
        """
    def evaluate(self, z: typing.SupportsInt, a: typing.SupportsInt, T9: typing.SupportsFloat) -> float:
        """
        Evaluate the composite partition function for given Z, A, and T9.
        """
    def evaluateDerivative(self, z: typing.SupportsInt, a: typing.SupportsInt, T9: typing.SupportsFloat) -> float:
        """
        Evaluate the derivative of the composite partition function for given Z, A, and T9.
        """
    def get_type(self) -> str:
        """
        Get the type of the partition function (should return 'Composite').
        """
    def supports(self, z: typing.SupportsInt, a: typing.SupportsInt) -> bool:
        """
        Check if the composite partition function supports given Z and A.
        """
class GroundStatePartitionFunction(PartitionFunction):
    def __init__(self) -> None:
        ...
    def evaluate(self, z: typing.SupportsInt, a: typing.SupportsInt, T9: typing.SupportsFloat) -> float:
        """
        Evaluate the ground state partition function for given Z, A, and T9.
        """
    def evaluateDerivative(self, z: typing.SupportsInt, a: typing.SupportsInt, T9: typing.SupportsFloat) -> float:
        """
        Evaluate the derivative of the ground state partition function for given Z, A, and T9.
        """
    def get_type(self) -> str:
        """
        Get the type of the partition function (should return 'GroundState').
        """
    def supports(self, z: typing.SupportsInt, a: typing.SupportsInt) -> bool:
        """
        Check if the ground state partition function supports given Z and A.
        """
class PartitionFunction:
    pass
class RauscherThielemannPartitionDataRecord:
    @property
    def a(self) -> int:
        """
        Mass number
        """
    @property
    def ground_state_spin(self) -> float:
        """
        Ground state spin
        """
    @property
    def normalized_g_values(self) -> float:
        """
        Normalized g-values for the first 24 energy levels
        """
    @property
    def z(self) -> int:
        """
        Atomic number
        """
class RauscherThielemannPartitionFunction(PartitionFunction):
    def __init__(self) -> None:
        ...
    def evaluate(self, z: typing.SupportsInt, a: typing.SupportsInt, T9: typing.SupportsFloat) -> float:
        """
        Evaluate the Rauscher-Thielemann partition function for given Z, A, and T9.
        """
    def evaluateDerivative(self, z: typing.SupportsInt, a: typing.SupportsInt, T9: typing.SupportsFloat) -> float:
        """
        Evaluate the derivative of the Rauscher-Thielemann partition function for given Z, A, and T9.
        """
    def get_type(self) -> str:
        """
        Get the type of the partition function (should return 'RauscherThielemann').
        """
    def supports(self, z: typing.SupportsInt, a: typing.SupportsInt) -> bool:
        """
        Check if the Rauscher-Thielemann partition function supports given Z and A.
        """
def basePartitionTypeToString(type: BasePartitionType) -> str:
    """
    Convert BasePartitionType to string.
    """
def stringToBasePartitionType(typeStr: str) -> BasePartitionType:
    """
    Convert string to BasePartitionType.
    """
GroundState: BasePartitionType  # value = <BasePartitionType.GroundState: 1>
RauscherThielemann: BasePartitionType  # value = <BasePartitionType.RauscherThielemann: 0>
