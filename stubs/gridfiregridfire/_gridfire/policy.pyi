"""
GridFire network policy bindings
"""
from __future__ import annotations
import collections.abc
import fourdst._phys.atomic
import fourdst._phys.composition
import gridfire._gridfire.engine
import gridfire._gridfire.engine.scratchpads
import gridfire._gridfire.partition
import gridfire._gridfire.reaction
import typing
__all__: list[str] = ['CNOChainPolicy', 'CNOIChainPolicy', 'CNOIIChainPolicy', 'CNOIIIChainPolicy', 'CNOIVChainPolicy', 'ConstructionResults', 'HotCNOChainPolicy', 'HotCNOIChainPolicy', 'HotCNOIIChainPolicy', 'HotCNOIIIChainPolicy', 'INITIALIZED_UNVERIFIED', 'INITIALIZED_VERIFIED', 'MISSING_KEY_REACTION', 'MISSING_KEY_SPECIES', 'MainSequencePolicy', 'MainSequenceReactionChainPolicy', 'MultiReactionChainPolicy', 'NetworkPolicy', 'NetworkPolicyStatus', 'ProtonProtonChainPolicy', 'ProtonProtonIChainPolicy', 'ProtonProtonIIChainPolicy', 'ProtonProtonIIIChainPolicy', 'ReactionChainPolicy', 'TemperatureDependentChainPolicy', 'TripleAlphaChainPolicy', 'UNINITIALIZED', 'network_policy_status_to_string']
class CNOChainPolicy(MultiReactionChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class CNOIChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class CNOIIChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class CNOIIIChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class CNOIVChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class ConstructionResults:
    @property
    def engine(self) -> gridfire._gridfire.engine.DynamicEngine:
        ...
    @property
    def scratch_blob(self) -> gridfire._gridfire.engine.scratchpads.StateBlob:
        ...
class HotCNOChainPolicy(MultiReactionChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class HotCNOIChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class HotCNOIIChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class HotCNOIIIChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class MainSequencePolicy(NetworkPolicy):
    @typing.overload
    def __init__(self, composition: fourdst._phys.composition.Composition) -> None:
        """
        Construct MainSequencePolicy from an existing composition.
        """
    @typing.overload
    def __init__(self, seed_species: collections.abc.Sequence[fourdst._phys.atomic.Species], mass_fractions: collections.abc.Sequence[typing.SupportsFloat]) -> None:
        """
        Construct MainSequencePolicy from seed species and mass fractions.
        """
    def construct(self) -> ConstructionResults:
        """
        Construct the network according to the policy.
        """
    def get_engine_stack(self) -> list[gridfire._gridfire.engine.DynamicEngine]:
        ...
    def get_engine_types_stack(self) -> list[gridfire._gridfire.engine.EngineTypes]:
        """
        Get the types of engines in the stack constructed by the network policy.
        """
    def get_partition_function(self) -> gridfire._gridfire.partition.PartitionFunction:
        ...
    def get_seed_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the set of seed reactions required by the network policy.
        """
    def get_seed_species(self) -> set[fourdst._phys.atomic.Species]:
        """
        Get the set of seed species required by the network policy.
        """
    def get_stack_scratch_blob(self) -> gridfire._gridfire.engine.scratchpads.StateBlob:
        ...
    def get_status(self) -> NetworkPolicyStatus:
        """
        Get the current status of the network policy.
        """
    def name(self) -> str:
        """
        Get the name of the network policy.
        """
class MainSequenceReactionChainPolicy(MultiReactionChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class MultiReactionChainPolicy(ReactionChainPolicy):
    pass
class NetworkPolicy:
    pass
class NetworkPolicyStatus:
    """
    Members:
    
      UNINITIALIZED
    
      INITIALIZED_UNVERIFIED
    
      MISSING_KEY_REACTION
    
      MISSING_KEY_SPECIES
    
      INITIALIZED_VERIFIED
    """
    INITIALIZED_UNVERIFIED: typing.ClassVar[NetworkPolicyStatus]  # value = <NetworkPolicyStatus.INITIALIZED_UNVERIFIED: 1>
    INITIALIZED_VERIFIED: typing.ClassVar[NetworkPolicyStatus]  # value = <NetworkPolicyStatus.INITIALIZED_VERIFIED: 4>
    MISSING_KEY_REACTION: typing.ClassVar[NetworkPolicyStatus]  # value = <NetworkPolicyStatus.MISSING_KEY_REACTION: 2>
    MISSING_KEY_SPECIES: typing.ClassVar[NetworkPolicyStatus]  # value = <NetworkPolicyStatus.MISSING_KEY_SPECIES: 3>
    UNINITIALIZED: typing.ClassVar[NetworkPolicyStatus]  # value = <NetworkPolicyStatus.UNINITIALIZED: 0>
    __members__: typing.ClassVar[dict[str, NetworkPolicyStatus]]  # value = {'UNINITIALIZED': <NetworkPolicyStatus.UNINITIALIZED: 0>, 'INITIALIZED_UNVERIFIED': <NetworkPolicyStatus.INITIALIZED_UNVERIFIED: 1>, 'MISSING_KEY_REACTION': <NetworkPolicyStatus.MISSING_KEY_REACTION: 2>, 'MISSING_KEY_SPECIES': <NetworkPolicyStatus.MISSING_KEY_SPECIES: 3>, 'INITIALIZED_VERIFIED': <NetworkPolicyStatus.INITIALIZED_VERIFIED: 4>}
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
class ProtonProtonChainPolicy(MultiReactionChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class ProtonProtonIChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class ProtonProtonIIChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class ProtonProtonIIIChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
class ReactionChainPolicy:
    pass
class TemperatureDependentChainPolicy(ReactionChainPolicy):
    pass
class TripleAlphaChainPolicy(TemperatureDependentChainPolicy):
    def __eq__(self, other: ReactionChainPolicy) -> bool:
        """
        Check equality with another ReactionChainPolicy.
        """
    def __hash__(self) -> int:
        ...
    def __init__(self) -> None:
        ...
    def __ne__(self, other: ReactionChainPolicy) -> bool:
        """
        Check inequality with another ReactionChainPolicy.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the reaction chain contains a reaction with the given ID.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the reaction chain contains the given reaction.
        """
    def get_reactions(self) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the ReactionSet representing this reaction chain.
        """
    def hash(self, seed: typing.SupportsInt) -> int:
        """
        Compute a hash value for the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
    @typing.overload
    def name(self) -> str:
        """
        Get the name of the reaction chain policy.
        """
def network_policy_status_to_string(status: NetworkPolicyStatus) -> str:
    """
    Convert a NetworkPolicyStatus enum value to its string representation.
    """
INITIALIZED_UNVERIFIED: NetworkPolicyStatus  # value = <NetworkPolicyStatus.INITIALIZED_UNVERIFIED: 1>
INITIALIZED_VERIFIED: NetworkPolicyStatus  # value = <NetworkPolicyStatus.INITIALIZED_VERIFIED: 4>
MISSING_KEY_REACTION: NetworkPolicyStatus  # value = <NetworkPolicyStatus.MISSING_KEY_REACTION: 2>
MISSING_KEY_SPECIES: NetworkPolicyStatus  # value = <NetworkPolicyStatus.MISSING_KEY_SPECIES: 3>
UNINITIALIZED: NetworkPolicyStatus  # value = <NetworkPolicyStatus.UNINITIALIZED: 0>
