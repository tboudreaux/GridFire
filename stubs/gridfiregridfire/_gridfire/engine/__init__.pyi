"""
Engine and Engine View bindings
"""
from __future__ import annotations
import collections.abc
import fourdst._phys.atomic
import fourdst._phys.composition
import gridfire._gridfire.io
import gridfire._gridfire.partition
import gridfire._gridfire.reaction
import gridfire._gridfire.screening
import gridfire._gridfire.type
import numpy
import numpy.typing
import typing
from . import diagnostics
from . import scratchpads
__all__: list[str] = ['ACTIVE', 'ADAPTIVE_ENGINE_VIEW', 'AdaptiveEngineView', 'BuildDepthType', 'DEFAULT', 'DEFINED_ENGINE_VIEW', 'DefinedEngineView', 'DynamicEngine', 'EQUILIBRIUM', 'Engine', 'EngineTypes', 'FILE_DEFINED_ENGINE_VIEW', 'FULL_SUCCESS', 'FifthOrder', 'FileDefinedEngineView', 'FourthOrder', 'Full', 'GRAPH_ENGINE', 'GraphEngine', 'INACTIVE_FLOW', 'MAX_ITERATIONS_REACHED', 'MULTISCALE_PARTITIONING_ENGINE_VIEW', 'MultiscalePartitioningEngineView', 'NONE', 'NOT_PRESENT', 'NO_SPECIES_TO_PRIME', 'NetworkBuildDepth', 'NetworkConstructionFlags', 'NetworkJacobian', 'NetworkPrimingEngineView', 'PRIMING_ENGINE_VIEW', 'PrimingReport', 'PrimingReportStatus', 'REACLIB', 'REACLIB_STRONG', 'REACLIB_WEAK', 'SecondOrder', 'Shallow', 'SparsityPattern', 'SpeciesStatus', 'StepDerivatives', 'ThirdOrder', 'WRL_BETA_MINUS', 'WRL_BETA_PLUS', 'WRL_ELECTRON_CAPTURE', 'WRL_POSITRON_CAPTURE', 'WRL_WEAK', 'build_nuclear_network', 'diagnostics', 'primeNetwork', 'regularize_jacobian', 'scratchpads']
class AdaptiveEngineView(DynamicEngine):
    def __init__(self, baseEngine: DynamicEngine) -> None:
        """
        Construct an adaptive engine view with a base engine.
        """
    def calculateEpsDerivatives(self, ctx: scratchpads.StateBlob, comp: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> ...:
        """
        Calculate deps/dT and deps/drho
        """
    def calculateMolarReactionFlow(self: DynamicEngine, ctx: scratchpads.StateBlob, reaction: ..., comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> float:
        """
        Calculate the molar reaction flow for a given reaction.
        """
    def calculateRHSAndEnergy(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> StepDerivatives:
        """
        Calculate the right-hand side (dY/dt) and energy generation rate.
        """
    def collectComposition(self, ctx: scratchpads.StateBlob, composition: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> fourdst._phys.composition.Composition:
        """
        Recursively collect composition from current engine and any sub engines if they exist.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> NetworkJacobian:
        """
        Generate the Jacobian matrix for the current state.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, activeSpecies: collections.abc.Sequence[fourdst._phys.atomic.Species]) -> NetworkJacobian:
        """
        Generate the jacobian matrix only for the subset of the matrix representing the active species.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, sparsityPattern: collections.abc.Sequence[tuple[typing.SupportsInt, typing.SupportsInt]]) -> NetworkJacobian:
        """
        Generate the jacobian matrix for the given sparsity pattern
        """
    def getBaseEngine(self) -> DynamicEngine:
        """
        Get the base engine associated with this adaptive engine view.
        """
    def getNetworkReactions(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the set of logical reactions in the network.
        """
    def getNetworkSpecies(self, arg0: scratchpads.StateBlob) -> list[fourdst._phys.atomic.Species]:
        """
        Get the list of species in the network.
        """
    def getScreeningModel(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.screening.ScreeningType:
        """
        Get the current screening model of the engine.
        """
    def getSpeciesDestructionTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the destruction timescales for each species in the network.
        """
    def getSpeciesIndex(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> int:
        """
        Get the index of a species in the network.
        """
    def getSpeciesStatus(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> SpeciesStatus:
        """
        Get the status of a species in the network.
        """
    def getSpeciesTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the timescales for each species in the network.
        """
    def primeEngine(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> PrimingReport:
        """
        Prime the engine with a NetIn object to prepare for calculations.
        """
    def project(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> fourdst._phys.composition.Composition:
        """
        Update the engine state based on the provided NetIn object.
        """
class BuildDepthType:
    pass
class DefinedEngineView(DynamicEngine):
    def __init__(self, peNames: collections.abc.Sequence[str], baseEngine: GraphEngine) -> None:
        """
        Construct a defined engine view with a list of tracked reactions and a base engine.
        """
    def calculateEpsDerivatives(self, ctx: scratchpads.StateBlob, comp: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> ...:
        """
        Calculate deps/dT and deps/drho
        """
    def calculateMolarReactionFlow(self: DynamicEngine, ctx: scratchpads.StateBlob, reaction: ..., comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> float:
        """
        Calculate the molar reaction flow for a given reaction.
        """
    def calculateRHSAndEnergy(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> StepDerivatives:
        """
        Calculate the right-hand side (dY/dt) and energy generation rate.
        """
    def collectComposition(self, ctx: scratchpads.StateBlob, composition: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> fourdst._phys.composition.Composition:
        """
        Recursively collect composition from current engine and any sub engines if they exist.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> NetworkJacobian:
        """
        Generate the Jacobian matrix for the current state.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, activeSpecies: collections.abc.Sequence[fourdst._phys.atomic.Species]) -> NetworkJacobian:
        """
        Generate the jacobian matrix only for the subset of the matrix representing the active species.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, sparsityPattern: collections.abc.Sequence[tuple[typing.SupportsInt, typing.SupportsInt]]) -> NetworkJacobian:
        """
        Generate the jacobian matrix for the given sparsity pattern
        """
    def getBaseEngine(self) -> DynamicEngine:
        """
        Get the base engine associated with this defined engine view.
        """
    def getNetworkReactions(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the set of logical reactions in the network.
        """
    def getNetworkSpecies(self, arg0: scratchpads.StateBlob) -> list[fourdst._phys.atomic.Species]:
        """
        Get the list of species in the network.
        """
    def getScreeningModel(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.screening.ScreeningType:
        """
        Get the current screening model of the engine.
        """
    def getSpeciesDestructionTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the destruction timescales for each species in the network.
        """
    def getSpeciesIndex(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> int:
        """
        Get the index of a species in the network.
        """
    def getSpeciesStatus(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> SpeciesStatus:
        """
        Get the status of a species in the network.
        """
    def getSpeciesTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the timescales for each species in the network.
        """
    def primeEngine(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> PrimingReport:
        """
        Prime the engine with a NetIn object to prepare for calculations.
        """
    def project(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> fourdst._phys.composition.Composition:
        """
        Update the engine state based on the provided NetIn object.
        """
class DynamicEngine:
    pass
class Engine:
    pass
class EngineTypes:
    """
    Members:
    
      GRAPH_ENGINE : The standard graph-based engine.
    
      ADAPTIVE_ENGINE_VIEW : An engine that adapts based on certain criteria.
    
      MULTISCALE_PARTITIONING_ENGINE_VIEW : An engine that partitions the system at multiple scales.
    
      PRIMING_ENGINE_VIEW : An engine that uses a priming strategy for simulations.
    
      DEFINED_ENGINE_VIEW : An engine defined by user specifications.
    
      FILE_DEFINED_ENGINE_VIEW : An engine defined through external files.
    """
    ADAPTIVE_ENGINE_VIEW: typing.ClassVar[EngineTypes]  # value = <EngineTypes.ADAPTIVE_ENGINE_VIEW: 1>
    DEFINED_ENGINE_VIEW: typing.ClassVar[EngineTypes]  # value = <EngineTypes.DEFINED_ENGINE_VIEW: 4>
    FILE_DEFINED_ENGINE_VIEW: typing.ClassVar[EngineTypes]  # value = <EngineTypes.FILE_DEFINED_ENGINE_VIEW: 5>
    GRAPH_ENGINE: typing.ClassVar[EngineTypes]  # value = <EngineTypes.GRAPH_ENGINE: 0>
    MULTISCALE_PARTITIONING_ENGINE_VIEW: typing.ClassVar[EngineTypes]  # value = <EngineTypes.MULTISCALE_PARTITIONING_ENGINE_VIEW: 2>
    PRIMING_ENGINE_VIEW: typing.ClassVar[EngineTypes]  # value = <EngineTypes.PRIMING_ENGINE_VIEW: 3>
    __members__: typing.ClassVar[dict[str, EngineTypes]]  # value = {'GRAPH_ENGINE': <EngineTypes.GRAPH_ENGINE: 0>, 'ADAPTIVE_ENGINE_VIEW': <EngineTypes.ADAPTIVE_ENGINE_VIEW: 1>, 'MULTISCALE_PARTITIONING_ENGINE_VIEW': <EngineTypes.MULTISCALE_PARTITIONING_ENGINE_VIEW: 2>, 'PRIMING_ENGINE_VIEW': <EngineTypes.PRIMING_ENGINE_VIEW: 3>, 'DEFINED_ENGINE_VIEW': <EngineTypes.DEFINED_ENGINE_VIEW: 4>, 'FILE_DEFINED_ENGINE_VIEW': <EngineTypes.FILE_DEFINED_ENGINE_VIEW: 5>}
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
    @typing.overload
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __repr__(self) -> str:
        """
        String representation of the EngineTypes.
        """
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
class FileDefinedEngineView(DefinedEngineView):
    def __init__(self, baseEngine: GraphEngine, fileName: str, parser: gridfire._gridfire.io.NetworkFileParser) -> None:
        """
        Construct a defined engine view from a file and a base engine.
        """
    def calculateEpsDerivatives(self, ctx: scratchpads.StateBlob, comp: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> ...:
        """
        Calculate deps/dT and deps/drho
        """
    def calculateMolarReactionFlow(self: DynamicEngine, ctx: scratchpads.StateBlob, reaction: ..., comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> float:
        """
        Calculate the molar reaction flow for a given reaction.
        """
    def calculateRHSAndEnergy(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> StepDerivatives:
        """
        Calculate the right-hand side (dY/dt) and energy generation rate.
        """
    def collectComposition(self, ctx: scratchpads.StateBlob, composition: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> fourdst._phys.composition.Composition:
        """
        Recursively collect composition from current engine and any sub engines if they exist.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> NetworkJacobian:
        """
        Generate the Jacobian matrix for the current state.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, activeSpecies: collections.abc.Sequence[fourdst._phys.atomic.Species]) -> NetworkJacobian:
        """
        Generate the jacobian matrix only for the subset of the matrix representing the active species.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, sparsityPattern: collections.abc.Sequence[tuple[typing.SupportsInt, typing.SupportsInt]]) -> NetworkJacobian:
        """
        Generate the jacobian matrix for the given sparsity pattern
        """
    def getBaseEngine(self) -> DynamicEngine:
        """
        Get the base engine associated with this file defined engine view.
        """
    def getNetworkFile(self) -> str:
        """
        Get the network file associated with this defined engine view.
        """
    def getNetworkReactions(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the set of logical reactions in the network.
        """
    def getNetworkSpecies(self, arg0: scratchpads.StateBlob) -> list[fourdst._phys.atomic.Species]:
        """
        Get the list of species in the network.
        """
    def getParser(self) -> gridfire._gridfire.io.NetworkFileParser:
        """
        Get the parser used for this defined engine view.
        """
    def getScreeningModel(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.screening.ScreeningType:
        """
        Get the current screening model of the engine.
        """
    def getSpeciesDestructionTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the destruction timescales for each species in the network.
        """
    def getSpeciesIndex(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> int:
        """
        Get the index of a species in the network.
        """
    def getSpeciesStatus(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> SpeciesStatus:
        """
        Get the status of a species in the network.
        """
    def getSpeciesTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the timescales for each species in the network.
        """
    def primeEngine(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> PrimingReport:
        """
        Prime the engine with a NetIn object to prepare for calculations.
        """
    def project(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> fourdst._phys.composition.Composition:
        """
        Update the engine state based on the provided NetIn object.
        """
class GraphEngine(DynamicEngine):
    @typing.overload
    def __init__(self, composition: fourdst._phys.composition.Composition, depth: gridfire._gridfire.engine.NetworkBuildDepth | typing.SupportsInt = ...) -> None:
        """
        Initialize GraphEngine with a composition and build depth.
        """
    @typing.overload
    def __init__(self, composition: fourdst._phys.composition.Composition, partitionFunction: gridfire._gridfire.partition.PartitionFunction, depth: gridfire._gridfire.engine.NetworkBuildDepth | typing.SupportsInt = ...) -> None:
        """
        Initialize GraphEngine with a composition, partition function and build depth.
        """
    @typing.overload
    def __init__(self, reactions: gridfire._gridfire.reaction.ReactionSet) -> None:
        """
        Initialize GraphEngine with a set of reactions.
        """
    def calculateEpsDerivatives(self, ctx: scratchpads.StateBlob, comp: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> ...:
        """
        Calculate deps/dT and deps/drho
        """
    def calculateMolarReactionFlow(self: DynamicEngine, ctx: scratchpads.StateBlob, reaction: ..., comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> float:
        """
        Calculate the molar reaction flow for a given reaction.
        """
    def calculateRHSAndEnergy(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> StepDerivatives:
        """
        Calculate the right-hand side (dY/dt) and energy generation rate.
        """
    def calculateReverseRate(self, reaction: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat, composition: ...) -> float:
        """
        Calculate the reverse rate for a given reaction at a specific temperature, density, and composition.
        """
    def calculateReverseRateTwoBody(self, reaction: ..., T9: typing.SupportsFloat, forwardRate: typing.SupportsFloat, expFactor: typing.SupportsFloat) -> float:
        """
        Calculate the reverse rate for a two-body reaction at a specific temperature.
        """
    def calculateReverseRateTwoBodyDerivative(self, reaction: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat, composition: fourdst._phys.composition.Composition, reverseRate: typing.SupportsFloat) -> float:
        """
        Calculate the derivative of the reverse rate for a two-body reaction at a specific temperature.
        """
    def collectComposition(self, ctx: scratchpads.StateBlob, composition: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> fourdst._phys.composition.Composition:
        """
        Recursively collect composition from current engine and any sub engines if they exist.
        """
    def exportToCSV(self, ctx: scratchpads.StateBlob, filename: str) -> None:
        """
        Export the network to a CSV file for analysis.
        """
    def exportToDot(self, ctx: scratchpads.StateBlob, filename: str) -> None:
        """
        Export the network to a DOT file for visualization.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> NetworkJacobian:
        """
        Generate the Jacobian matrix for the current state.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, activeSpecies: collections.abc.Sequence[fourdst._phys.atomic.Species]) -> NetworkJacobian:
        """
        Generate the jacobian matrix only for the subset of the matrix representing the active species.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, sparsityPattern: collections.abc.Sequence[tuple[typing.SupportsInt, typing.SupportsInt]]) -> NetworkJacobian:
        """
        Generate the jacobian matrix for the given sparsity pattern
        """
    def getNetworkReactions(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the set of logical reactions in the network.
        """
    def getNetworkSpecies(self, arg0: scratchpads.StateBlob) -> list[fourdst._phys.atomic.Species]:
        """
        Get the list of species in the network.
        """
    def getPartitionFunction(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.partition.PartitionFunction:
        """
        Get the partition function used by the engine.
        """
    def getScreeningModel(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.screening.ScreeningType:
        """
        Get the current screening model of the engine.
        """
    @typing.overload
    def getSpeciesDestructionTimescales(self, ctx: scratchpads.StateBlob, composition: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, activeReactions: gridfire._gridfire.reaction.ReactionSet) -> ...:
        ...
    @typing.overload
    def getSpeciesDestructionTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the destruction timescales for each species in the network.
        """
    def getSpeciesIndex(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> int:
        """
        Get the index of a species in the network.
        """
    def getSpeciesStatus(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> SpeciesStatus:
        """
        Get the status of a species in the network.
        """
    @typing.overload
    def getSpeciesTimescales(self, ctx: scratchpads.StateBlob, composition: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, activeReactions: gridfire._gridfire.reaction.ReactionSet) -> ...:
        ...
    @typing.overload
    def getSpeciesTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the timescales for each species in the network.
        """
    def involvesSpecies(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> bool:
        """
        Check if a given species is involved in the network.
        """
    def isPrecomputationEnabled(self, arg0: scratchpads.StateBlob) -> bool:
        """
        Check if precomputation is enabled for the engine.
        """
    def isUsingReverseReactions(self, arg0: scratchpads.StateBlob) -> bool:
        """
        Check if the engine is using reverse reactions.
        """
    def primeEngine(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> PrimingReport:
        """
        Prime the engine with a NetIn object to prepare for calculations.
        """
    def project(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> fourdst._phys.composition.Composition:
        """
        Update the engine state based on the provided NetIn object.
        """
class MultiscalePartitioningEngineView(DynamicEngine):
    def __init__(self, baseEngine: GraphEngine) -> None:
        """
        Construct a multiscale partitioning engine view with a base engine.
        """
    def calculateEpsDerivatives(self, ctx: scratchpads.StateBlob, comp: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> ...:
        """
        Calculate deps/dT and deps/drho
        """
    def calculateMolarReactionFlow(self: DynamicEngine, ctx: scratchpads.StateBlob, reaction: ..., comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> float:
        """
        Calculate the molar reaction flow for a given reaction.
        """
    def calculateRHSAndEnergy(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> StepDerivatives:
        """
        Calculate the right-hand side (dY/dt) and energy generation rate.
        """
    def collectComposition(self, ctx: scratchpads.StateBlob, composition: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> fourdst._phys.composition.Composition:
        """
        Recursively collect composition from current engine and any sub engines if they exist.
        """
    def exportToDot(self, ctx: scratchpads.StateBlob, filename: str, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> None:
        """
        Export the network to a DOT file for visualization.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> NetworkJacobian:
        """
        Generate the Jacobian matrix for the current state.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, activeSpecies: collections.abc.Sequence[fourdst._phys.atomic.Species]) -> NetworkJacobian:
        """
        Generate the jacobian matrix only for the subset of the matrix representing the active species.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, sparsityPattern: collections.abc.Sequence[tuple[typing.SupportsInt, typing.SupportsInt]]) -> NetworkJacobian:
        """
        Generate the jacobian matrix for the given sparsity pattern
        """
    def getBaseEngine(self) -> DynamicEngine:
        """
        Get the base engine associated with this multiscale partitioning engine view.
        """
    def getDynamicSpecies(self: scratchpads.StateBlob) -> list[fourdst._phys.atomic.Species]:
        """
        Get the list of dynamic species in the network.
        """
    def getFastSpecies(self, arg0: scratchpads.StateBlob) -> list[fourdst._phys.atomic.Species]:
        """
        Get the list of fast species in the network.
        """
    def getNetworkReactions(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the set of logical reactions in the network.
        """
    def getNetworkSpecies(self, arg0: scratchpads.StateBlob) -> list[fourdst._phys.atomic.Species]:
        """
        Get the list of species in the network.
        """
    def getNormalizedEquilibratedComposition(self, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> fourdst._phys.composition.Composition:
        """
        Get the normalized equilibrated composition for the algebraic species.
        """
    def getScreeningModel(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.screening.ScreeningType:
        """
        Get the current screening model of the engine.
        """
    def getSpeciesDestructionTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the destruction timescales for each species in the network.
        """
    def getSpeciesIndex(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> int:
        """
        Get the index of a species in the network.
        """
    def getSpeciesStatus(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> SpeciesStatus:
        """
        Get the status of a species in the network.
        """
    def getSpeciesTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the timescales for each species in the network.
        """
    def involvesSpecies(self: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> bool:
        """
        Check if a given species is involved in the network (in either the algebraic or dynamic set).
        """
    def involvesSpeciesInDynamic(self: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> bool:
        """
        Check if a given species is involved in the network's dynamic set.
        """
    def involvesSpeciesInQSE(self: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> bool:
        """
        Check if a given species is involved in the network's algebraic set.
        """
    def partitionNetwork(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> fourdst._phys.composition.Composition:
        """
        Partition the network based on species timescales and connectivity.
        """
    def primeEngine(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> PrimingReport:
        """
        Prime the engine with a NetIn object to prepare for calculations.
        """
    def project(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> fourdst._phys.composition.Composition:
        """
        Update the engine state based on the provided NetIn object.
        """
class NetworkBuildDepth:
    """
    Members:
    
      Full : Full network build depth
    
      Shallow : Shallow network build depth
    
      SecondOrder : Second order network build depth
    
      ThirdOrder : Third order network build depth
    
      FourthOrder : Fourth order network build depth
    
      FifthOrder : Fifth order network build depth
    """
    FifthOrder: typing.ClassVar[NetworkBuildDepth]  # value = <NetworkBuildDepth.FifthOrder: 5>
    FourthOrder: typing.ClassVar[NetworkBuildDepth]  # value = <NetworkBuildDepth.FourthOrder: 4>
    Full: typing.ClassVar[NetworkBuildDepth]  # value = <NetworkBuildDepth.Full: -1>
    SecondOrder: typing.ClassVar[NetworkBuildDepth]  # value = <NetworkBuildDepth.SecondOrder: 2>
    Shallow: typing.ClassVar[NetworkBuildDepth]  # value = <NetworkBuildDepth.Shallow: 1>
    ThirdOrder: typing.ClassVar[NetworkBuildDepth]  # value = <NetworkBuildDepth.ThirdOrder: 3>
    __members__: typing.ClassVar[dict[str, NetworkBuildDepth]]  # value = {'Full': <NetworkBuildDepth.Full: -1>, 'Shallow': <NetworkBuildDepth.Shallow: 1>, 'SecondOrder': <NetworkBuildDepth.SecondOrder: 2>, 'ThirdOrder': <NetworkBuildDepth.ThirdOrder: 3>, 'FourthOrder': <NetworkBuildDepth.FourthOrder: 4>, 'FifthOrder': <NetworkBuildDepth.FifthOrder: 5>}
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
class NetworkConstructionFlags:
    """
    Members:
    
      NONE : No special construction flags.
    
      REACLIB_STRONG : Include strong reactions from reaclib.
    
      WRL_BETA_MINUS : Include beta-minus decay reactions from weak rate library.
    
      WRL_BETA_PLUS : Include beta-plus decay reactions from weak rate library.
    
      WRL_ELECTRON_CAPTURE : Include electron capture reactions from weak rate library.
    
      WRL_POSITRON_CAPTURE : Include positron capture reactions from weak rate library.
    
      REACLIB_WEAK : Include weak reactions from reaclib.
    
      WRL_WEAK : Include all weak reactions from weak rate library.
    
      REACLIB : Include all reactions from reaclib.
    
      DEFAULT : Default construction flags (Reaclib strong and weak).
    """
    DEFAULT: typing.ClassVar[NetworkConstructionFlags]  # value = <NetworkConstructionFlags.REACLIB: 33>
    NONE: typing.ClassVar[NetworkConstructionFlags]  # value = <NetworkConstructionFlags.NONE: 0>
    REACLIB: typing.ClassVar[NetworkConstructionFlags]  # value = <NetworkConstructionFlags.REACLIB: 33>
    REACLIB_STRONG: typing.ClassVar[NetworkConstructionFlags]  # value = <NetworkConstructionFlags.REACLIB_STRONG: 1>
    REACLIB_WEAK: typing.ClassVar[NetworkConstructionFlags]  # value = <NetworkConstructionFlags.REACLIB_WEAK: 32>
    WRL_BETA_MINUS: typing.ClassVar[NetworkConstructionFlags]  # value = <NetworkConstructionFlags.WRL_BETA_MINUS: 2>
    WRL_BETA_PLUS: typing.ClassVar[NetworkConstructionFlags]  # value = <NetworkConstructionFlags.WRL_BETA_PLUS: 4>
    WRL_ELECTRON_CAPTURE: typing.ClassVar[NetworkConstructionFlags]  # value = <NetworkConstructionFlags.WRL_ELECTRON_CAPTURE: 8>
    WRL_POSITRON_CAPTURE: typing.ClassVar[NetworkConstructionFlags]  # value = <NetworkConstructionFlags.WRL_POSITRON_CAPTURE: 16>
    WRL_WEAK: typing.ClassVar[NetworkConstructionFlags]  # value = <NetworkConstructionFlags.WRL_WEAK: 30>
    __members__: typing.ClassVar[dict[str, NetworkConstructionFlags]]  # value = {'NONE': <NetworkConstructionFlags.NONE: 0>, 'REACLIB_STRONG': <NetworkConstructionFlags.REACLIB_STRONG: 1>, 'WRL_BETA_MINUS': <NetworkConstructionFlags.WRL_BETA_MINUS: 2>, 'WRL_BETA_PLUS': <NetworkConstructionFlags.WRL_BETA_PLUS: 4>, 'WRL_ELECTRON_CAPTURE': <NetworkConstructionFlags.WRL_ELECTRON_CAPTURE: 8>, 'WRL_POSITRON_CAPTURE': <NetworkConstructionFlags.WRL_POSITRON_CAPTURE: 16>, 'REACLIB_WEAK': <NetworkConstructionFlags.REACLIB_WEAK: 32>, 'WRL_WEAK': <NetworkConstructionFlags.WRL_WEAK: 30>, 'REACLIB': <NetworkConstructionFlags.REACLIB: 33>, 'DEFAULT': <NetworkConstructionFlags.REACLIB: 33>}
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
    @typing.overload
    def __repr__(self) -> str:
        ...
    @typing.overload
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
class NetworkJacobian:
    @typing.overload
    def __getitem__(self, key: tuple[fourdst._phys.atomic.Species, fourdst._phys.atomic.Species]) -> float:
        """
        Get an entry from the Jacobian matrix using species identifiers.
        """
    @typing.overload
    def __getitem__(self, key: tuple[typing.SupportsInt, typing.SupportsInt]) -> float:
        """
        Get an entry from the Jacobian matrix using indices.
        """
    @typing.overload
    def __setitem__(self, key: tuple[fourdst._phys.atomic.Species, fourdst._phys.atomic.Species], value: typing.SupportsFloat) -> None:
        """
        Set an entry in the Jacobian matrix using species identifiers.
        """
    @typing.overload
    def __setitem__(self, key: tuple[typing.SupportsInt, typing.SupportsInt], value: typing.SupportsFloat) -> None:
        """
        Set an entry in the Jacobian matrix using indices.
        """
    def data(self) -> ...:
        """
        Get the underlying sparse matrix data.
        """
    def infs(self) -> list[tuple[tuple[fourdst._phys.atomic.Species, fourdst._phys.atomic.Species], float]]:
        """
        Get all infinite entries in the Jacobian matrix.
        """
    def mapping(self) -> dict[fourdst._phys.atomic.Species, int]:
        """
        Get the species-to-index mapping.
        """
    def nans(self) -> list[tuple[tuple[fourdst._phys.atomic.Species, fourdst._phys.atomic.Species], float]]:
        """
        Get all NaN entries in the Jacobian matrix.
        """
    def nnz(self) -> int:
        """
        Get the number of non-zero entries in the Jacobian matrix.
        """
    def rank(self) -> int:
        """
        Get the rank of the Jacobian matrix.
        """
    def shape(self) -> tuple[int, int]:
        """
        Get the shape of the Jacobian matrix as (rows, columns).
        """
    def singular(self) -> bool:
        """
        Check if the Jacobian matrix is singular.
        """
    def to_csv(self, filename: str) -> None:
        """
        Export the Jacobian matrix to a CSV file.
        """
    def to_numpy(self) -> numpy.typing.NDArray[numpy.float64]:
        """
        Convert the Jacobian matrix to a NumPy array.
        """
class NetworkPrimingEngineView(DefinedEngineView):
    @typing.overload
    def __init__(self, ctx: scratchpads.StateBlob, primingSymbol: str, baseEngine: GraphEngine) -> None:
        """
        Construct a priming engine view with a priming symbol and a base engine.
        """
    @typing.overload
    def __init__(self, ctx: scratchpads.StateBlob, primingSpecies: fourdst._phys.atomic.Species, baseEngine: GraphEngine) -> None:
        """
        Construct a priming engine view with a priming species and a base engine.
        """
    def calculateEpsDerivatives(self, ctx: scratchpads.StateBlob, comp: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> ...:
        """
        Calculate deps/dT and deps/drho
        """
    def calculateMolarReactionFlow(self: DynamicEngine, ctx: scratchpads.StateBlob, reaction: ..., comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> float:
        """
        Calculate the molar reaction flow for a given reaction.
        """
    def calculateRHSAndEnergy(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> StepDerivatives:
        """
        Calculate the right-hand side (dY/dt) and energy generation rate.
        """
    def collectComposition(self, ctx: scratchpads.StateBlob, composition: ..., T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> fourdst._phys.composition.Composition:
        """
        Recursively collect composition from current engine and any sub engines if they exist.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> NetworkJacobian:
        """
        Generate the Jacobian matrix for the current state.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, activeSpecies: collections.abc.Sequence[fourdst._phys.atomic.Species]) -> NetworkJacobian:
        """
        Generate the jacobian matrix only for the subset of the matrix representing the active species.
        """
    @typing.overload
    def generateJacobianMatrix(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat, sparsityPattern: collections.abc.Sequence[tuple[typing.SupportsInt, typing.SupportsInt]]) -> NetworkJacobian:
        """
        Generate the jacobian matrix for the given sparsity pattern
        """
    def getBaseEngine(self) -> DynamicEngine:
        """
        Get the base engine associated with this priming engine view.
        """
    def getNetworkReactions(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.reaction.ReactionSet:
        """
        Get the set of logical reactions in the network.
        """
    def getNetworkSpecies(self, arg0: scratchpads.StateBlob) -> list[fourdst._phys.atomic.Species]:
        """
        Get the list of species in the network.
        """
    def getScreeningModel(self, arg0: scratchpads.StateBlob) -> gridfire._gridfire.screening.ScreeningType:
        """
        Get the current screening model of the engine.
        """
    def getSpeciesDestructionTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the destruction timescales for each species in the network.
        """
    def getSpeciesIndex(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> int:
        """
        Get the index of a species in the network.
        """
    def getSpeciesStatus(self, ctx: scratchpads.StateBlob, species: fourdst._phys.atomic.Species) -> SpeciesStatus:
        """
        Get the status of a species in the network.
        """
    def getSpeciesTimescales(self: DynamicEngine, ctx: scratchpads.StateBlob, comp: fourdst._phys.composition.Composition, T9: typing.SupportsFloat, rho: typing.SupportsFloat) -> dict[fourdst._phys.atomic.Species, float]:
        """
        Get the timescales for each species in the network.
        """
    def primeEngine(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> PrimingReport:
        """
        Prime the engine with a NetIn object to prepare for calculations.
        """
    def project(self, ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn) -> fourdst._phys.composition.Composition:
        """
        Update the engine state based on the provided NetIn object.
        """
class PrimingReport:
    def __repr__(self) -> str:
        ...
    @property
    def primedComposition(self) -> fourdst._phys.composition.Composition:
        """
        The composition after priming.
        """
    @property
    def status(self) -> PrimingReportStatus:
        """
        Status message from the priming process.
        """
    @property
    def success(self) -> bool:
        """
        Indicates if the priming was successful.
        """
class PrimingReportStatus:
    """
    Members:
    
      FULL_SUCCESS : Priming was full successful.
    
      NO_SPECIES_TO_PRIME : Solver Failed to converge during priming.
    
      MAX_ITERATIONS_REACHED : Engine has already been primed.
    """
    FULL_SUCCESS: typing.ClassVar[PrimingReportStatus]  # value = <PrimingReportStatus.FULL_SUCCESS: 0>
    MAX_ITERATIONS_REACHED: typing.ClassVar[PrimingReportStatus]  # value = <PrimingReportStatus.MAX_ITERATIONS_REACHED: 1>
    NO_SPECIES_TO_PRIME: typing.ClassVar[PrimingReportStatus]  # value = <PrimingReportStatus.NO_SPECIES_TO_PRIME: 2>
    __members__: typing.ClassVar[dict[str, PrimingReportStatus]]  # value = {'FULL_SUCCESS': <PrimingReportStatus.FULL_SUCCESS: 0>, 'NO_SPECIES_TO_PRIME': <PrimingReportStatus.NO_SPECIES_TO_PRIME: 2>, 'MAX_ITERATIONS_REACHED': <PrimingReportStatus.MAX_ITERATIONS_REACHED: 1>}
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
    @typing.overload
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __repr__(self) -> str:
        """
        String representation of the PrimingReport.
        """
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
class SparsityPattern:
    pass
class SpeciesStatus:
    """
    Members:
    
      ACTIVE : Species is active in the network.
    
      EQUILIBRIUM : Species is in equilibrium.
    
      INACTIVE_FLOW : Species is inactive due to flow.
    
      NOT_PRESENT : Species is not present in the network.
    """
    ACTIVE: typing.ClassVar[SpeciesStatus]  # value = <SpeciesStatus.ACTIVE: 0>
    EQUILIBRIUM: typing.ClassVar[SpeciesStatus]  # value = <SpeciesStatus.EQUILIBRIUM: 1>
    INACTIVE_FLOW: typing.ClassVar[SpeciesStatus]  # value = <SpeciesStatus.INACTIVE_FLOW: 2>
    NOT_PRESENT: typing.ClassVar[SpeciesStatus]  # value = <SpeciesStatus.NOT_PRESENT: 3>
    __members__: typing.ClassVar[dict[str, SpeciesStatus]]  # value = {'ACTIVE': <SpeciesStatus.ACTIVE: 0>, 'EQUILIBRIUM': <SpeciesStatus.EQUILIBRIUM: 1>, 'INACTIVE_FLOW': <SpeciesStatus.INACTIVE_FLOW: 2>, 'NOT_PRESENT': <SpeciesStatus.NOT_PRESENT: 3>}
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
    @typing.overload
    def __repr__(self) -> str:
        ...
    @typing.overload
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
class StepDerivatives:
    @property
    def dYdt(self) -> dict[fourdst._phys.atomic.Species, float]:
        """
        The right-hand side (dY/dt) of the ODE system.
        """
    @property
    def energy(self) -> float:
        """
        The energy generation rate.
        """
def build_nuclear_network(composition: ..., weakInterpolator: ..., maxLayers: gridfire._gridfire.engine.NetworkBuildDepth | typing.SupportsInt = ..., ReactionTypes: NetworkConstructionFlags = ...) -> gridfire._gridfire.reaction.ReactionSet:
    """
    Build a nuclear network from a composition using all archived reaction data.
    """
def primeNetwork(ctx: scratchpads.StateBlob, netIn: gridfire._gridfire.type.NetIn, engine: ..., ignoredReactionTypes: collections.abc.Sequence[...] | None = None) -> PrimingReport:
    """
    Prime a network with a short timescale ignition
    """
def regularize_jacobian(jacobian: NetworkJacobian, composition: fourdst._phys.composition.Composition) -> NetworkJacobian:
    """
    regularize_jacobian
    """
ACTIVE: SpeciesStatus  # value = <SpeciesStatus.ACTIVE: 0>
ADAPTIVE_ENGINE_VIEW: EngineTypes  # value = <EngineTypes.ADAPTIVE_ENGINE_VIEW: 1>
DEFAULT: NetworkConstructionFlags  # value = <NetworkConstructionFlags.REACLIB: 33>
DEFINED_ENGINE_VIEW: EngineTypes  # value = <EngineTypes.DEFINED_ENGINE_VIEW: 4>
EQUILIBRIUM: SpeciesStatus  # value = <SpeciesStatus.EQUILIBRIUM: 1>
FILE_DEFINED_ENGINE_VIEW: EngineTypes  # value = <EngineTypes.FILE_DEFINED_ENGINE_VIEW: 5>
FULL_SUCCESS: PrimingReportStatus  # value = <PrimingReportStatus.FULL_SUCCESS: 0>
FifthOrder: NetworkBuildDepth  # value = <NetworkBuildDepth.FifthOrder: 5>
FourthOrder: NetworkBuildDepth  # value = <NetworkBuildDepth.FourthOrder: 4>
Full: NetworkBuildDepth  # value = <NetworkBuildDepth.Full: -1>
GRAPH_ENGINE: EngineTypes  # value = <EngineTypes.GRAPH_ENGINE: 0>
INACTIVE_FLOW: SpeciesStatus  # value = <SpeciesStatus.INACTIVE_FLOW: 2>
MAX_ITERATIONS_REACHED: PrimingReportStatus  # value = <PrimingReportStatus.MAX_ITERATIONS_REACHED: 1>
MULTISCALE_PARTITIONING_ENGINE_VIEW: EngineTypes  # value = <EngineTypes.MULTISCALE_PARTITIONING_ENGINE_VIEW: 2>
NONE: NetworkConstructionFlags  # value = <NetworkConstructionFlags.NONE: 0>
NOT_PRESENT: SpeciesStatus  # value = <SpeciesStatus.NOT_PRESENT: 3>
NO_SPECIES_TO_PRIME: PrimingReportStatus  # value = <PrimingReportStatus.NO_SPECIES_TO_PRIME: 2>
PRIMING_ENGINE_VIEW: EngineTypes  # value = <EngineTypes.PRIMING_ENGINE_VIEW: 3>
REACLIB: NetworkConstructionFlags  # value = <NetworkConstructionFlags.REACLIB: 33>
REACLIB_STRONG: NetworkConstructionFlags  # value = <NetworkConstructionFlags.REACLIB_STRONG: 1>
REACLIB_WEAK: NetworkConstructionFlags  # value = <NetworkConstructionFlags.REACLIB_WEAK: 32>
SecondOrder: NetworkBuildDepth  # value = <NetworkBuildDepth.SecondOrder: 2>
Shallow: NetworkBuildDepth  # value = <NetworkBuildDepth.Shallow: 1>
ThirdOrder: NetworkBuildDepth  # value = <NetworkBuildDepth.ThirdOrder: 3>
WRL_BETA_MINUS: NetworkConstructionFlags  # value = <NetworkConstructionFlags.WRL_BETA_MINUS: 2>
WRL_BETA_PLUS: NetworkConstructionFlags  # value = <NetworkConstructionFlags.WRL_BETA_PLUS: 4>
WRL_ELECTRON_CAPTURE: NetworkConstructionFlags  # value = <NetworkConstructionFlags.WRL_ELECTRON_CAPTURE: 8>
WRL_POSITRON_CAPTURE: NetworkConstructionFlags  # value = <NetworkConstructionFlags.WRL_POSITRON_CAPTURE: 16>
WRL_WEAK: NetworkConstructionFlags  # value = <NetworkConstructionFlags.WRL_WEAK: 30>
