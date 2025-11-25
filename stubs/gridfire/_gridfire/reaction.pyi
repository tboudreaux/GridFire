"""
GridFire reaction bindings
"""
from __future__ import annotations
import collections.abc
import fourdst._phys.atomic
import fourdst._phys.composition
import typing
__all__: list[str] = ['LogicalReaclibReaction', 'RateCoefficientSet', 'ReaclibReaction', 'ReactionSet', 'get_all_reactions', 'packReactionSet']
class LogicalReaclibReaction(ReaclibReaction):
    @typing.overload
    def __init__(self, reactions: collections.abc.Sequence[ReaclibReaction]) -> None:
        """
        Construct a LogicalReaclibReaction from a vector of ReaclibReaction objects.
        """
    @typing.overload
    def __init__(self, reactions: collections.abc.Sequence[ReaclibReaction], is_reverse: bool) -> None:
        """
        Construct a LogicalReaclibReaction from a vector of ReaclibReaction objects.
        """
    def __len__(self) -> int:
        """
        Overload len() to return the number of source rates.
        """
    def add_reaction(self, reaction: ReaclibReaction) -> None:
        """
        Add another Reaction source to this logical reaction.
        """
    def calculate_forward_rate_log_derivative(self, T9: typing.SupportsFloat, rho: typing.SupportsFloat, Ye: typing.SupportsFloat, mue: typing.SupportsFloat, Composition: fourdst._phys.composition.Composition) -> float:
        """
        Calculate the forward rate log derivative at a given temperature T9 (in units of 10^9 K).
        """
    def calculate_rate(self, T9: typing.SupportsFloat, rho: typing.SupportsFloat, Ye: typing.SupportsFloat, mue: typing.SupportsFloat, Y: collections.abc.Sequence[typing.SupportsFloat], index_to_species_map: collections.abc.Mapping[typing.SupportsInt, fourdst._phys.atomic.Species]) -> float:
        """
        Calculate the reaction rate at a given temperature T9 (in units of 10^9 K). Note that for a reaclib reaction only T9 is actually used, all other parameters are there for interface compatibility.
        """
    def size(self) -> int:
        """
        Get the number of source rates contributing to this logical reaction.
        """
    def sources(self) -> list[str]:
        """
        Get the list of source labels for the aggregated rates.
        """
class RateCoefficientSet:
    def __init__(self, a0: typing.SupportsFloat, a1: typing.SupportsFloat, a2: typing.SupportsFloat, a3: typing.SupportsFloat, a4: typing.SupportsFloat, a5: typing.SupportsFloat, a6: typing.SupportsFloat) -> None:
        """
        Construct a RateCoefficientSet with the given parameters.
        """
class ReaclibReaction:
    __hash__: typing.ClassVar[None] = None
    def __eq__(self, arg0: ReaclibReaction) -> bool:
        """
        Equality operator for reactions based on their IDs.
        """
    def __init__(self, id: str, peName: str, chapter: typing.SupportsInt, reactants: collections.abc.Sequence[fourdst._phys.atomic.Species], products: collections.abc.Sequence[fourdst._phys.atomic.Species], qValue: typing.SupportsFloat, label: str, sets: RateCoefficientSet, reverse: bool = False) -> None:
        """
        Construct a Reaction with the given parameters.
        """
    def __neq__(self, arg0: ReaclibReaction) -> bool:
        """
        Inequality operator for reactions based on their IDs.
        """
    def __repr__(self) -> str:
        ...
    def all_species(self) -> set[fourdst._phys.atomic.Species]:
        """
        Get all species involved in the reaction (both reactants and products) as a set.
        """
    def calculate_rate(self, T9: typing.SupportsFloat, rho: typing.SupportsFloat, Y: collections.abc.Sequence[typing.SupportsFloat]) -> float:
        """
        Calculate the reaction rate at a given temperature T9 (in units of 10^9 K).
        """
    def chapter(self) -> int:
        """
        Get the REACLIB chapter number defining the reaction structure.
        """
    def contains(self, species: fourdst._phys.atomic.Species) -> bool:
        """
        Check if the reaction contains a specific species.
        """
    def contains_product(self, arg0: fourdst._phys.atomic.Species) -> bool:
        """
        Check if the reaction contains a specific product species.
        """
    def contains_reactant(self, arg0: fourdst._phys.atomic.Species) -> bool:
        """
        Check if the reaction contains a specific reactant species.
        """
    def excess_energy(self) -> float:
        """
        Calculate the excess energy from the mass difference of reactants and products.
        """
    def hash(self, seed: typing.SupportsInt = 0) -> int:
        """
        Compute a hash for the reaction based on its ID.
        """
    def id(self) -> str:
        """
        Get the unique identifier of the reaction.
        """
    def is_reverse(self) -> bool:
        """
        Check if this is a reverse reaction rate.
        """
    def num_species(self) -> int:
        """
        Count the number of species in the reaction.
        """
    def peName(self) -> str:
        """
        Get the reaction name in (projectile, ejectile) notation (e.g., 'p(p,g)d').
        """
    def product_species(self) -> set[fourdst._phys.atomic.Species]:
        """
        Get the product species of the reaction as a set.
        """
    def products(self) -> list[fourdst._phys.atomic.Species]:
        """
        Get a list of product species in the reaction.
        """
    def qValue(self) -> float:
        """
        Get the Q-value of the reaction in MeV.
        """
    def rateCoefficients(self) -> RateCoefficientSet:
        """
        get the set of rate coefficients.
        """
    def reactant_species(self) -> set[fourdst._phys.atomic.Species]:
        """
        Get the reactant species of the reaction as a set.
        """
    def reactants(self) -> list[fourdst._phys.atomic.Species]:
        """
        Get a list of reactant species in the reaction.
        """
    def sourceLabel(self) -> str:
        """
        Get the source label for the rate data (e.g., 'wc12w', 'st08').
        """
    @typing.overload
    def stoichiometry(self, species: fourdst._phys.atomic.Species) -> int:
        """
        Get the stoichiometry of the reaction as a map from species to their coefficients.
        """
    @typing.overload
    def stoichiometry(self) -> dict[fourdst._phys.atomic.Species, int]:
        """
        Get the stoichiometry of the reaction as a map from species to their coefficients.
        """
class ReactionSet:
    __hash__: typing.ClassVar[None] = None
    @staticmethod
    def from_clones(reactions: collections.abc.Sequence[...]) -> ReactionSet:
        """
        Create a ReactionSet that takes ownership of the reactions by cloning the input reactions.
        """
    def __eq__(self, LogicalReactionSet: ReactionSet) -> bool:
        """
        Equality operator for LogicalReactionSets based on their contents.
        """
    def __getitem__(self, index: typing.SupportsInt) -> ...:
        """
        Get a LogicalReaclibReaction by index.
        """
    def __getitem___(self, id: str) -> ...:
        """
        Get a LogicalReaclibReaction by its ID.
        """
    @typing.overload
    def __init__(self, reactions: collections.abc.Sequence[...]) -> None:
        """
        Construct a LogicalReactionSet from a vector of LogicalReaclibReaction objects.
        """
    @typing.overload
    def __init__(self) -> None:
        """
        Default constructor for an empty LogicalReactionSet.
        """
    @typing.overload
    def __init__(self, other: ReactionSet) -> None:
        """
        Copy constructor for LogicalReactionSet.
        """
    def __len__(self) -> int:
        """
        Overload len() to return the number of LogicalReactions.
        """
    def __ne__(self, LogicalReactionSet: ReactionSet) -> bool:
        """
        Inequality operator for LogicalReactionSets based on their contents.
        """
    def __repr__(self) -> str:
        ...
    def add_reaction(self, reaction: ...) -> None:
        """
        Add a LogicalReaclibReaction to the set.
        """
    def clear(self) -> None:
        """
        Remove all LogicalReactions from the set.
        """
    @typing.overload
    def contains(self, id: str) -> bool:
        """
        Check if the set contains a specific LogicalReaclibReaction.
        """
    @typing.overload
    def contains(self, reaction: ...) -> bool:
        """
        Check if the set contains a specific Reaction.
        """
    def contains_product(self, species: fourdst._phys.atomic.Species) -> bool:
        """
        Check if any reaction in the set has the species as a product.
        """
    def contains_reactant(self, species: fourdst._phys.atomic.Species) -> bool:
        """
        Check if any reaction in the set has the species as a reactant.
        """
    def contains_species(self, species: fourdst._phys.atomic.Species) -> bool:
        """
        Check if any reaction in the set involves the given species.
        """
    def getReactionSetSpecies(self) -> set[fourdst._phys.atomic.Species]:
        """
        Get all species involved in the reactions of the set as a set of Species objects.
        """
    def hash(self, seed: typing.SupportsInt = 0) -> int:
        """
        Compute a hash for the LogicalReactionSet based on its contents.
        """
    def remove_reaction(self, reaction: ...) -> None:
        """
        Remove a LogicalReaclibReaction from the set.
        """
    def size(self) -> int:
        """
        Get the number of LogicalReactions in the set.
        """
def get_all_reactions() -> ReactionSet:
    """
    Get all reactions from the REACLIB database.
    """
def packReactionSet(reactionSet: ReactionSet) -> ReactionSet:
    """
    Convert a ReactionSet to a LogicalReactionSet by aggregating reactions with the same peName.
    """
