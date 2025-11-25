"""
GridFire exceptions bindings
"""
from __future__ import annotations
__all__: list[str] = ['BadCollectionError', 'BadRHSEngineError', 'CVODESolverFailureError', 'DebugException', 'EngineError', 'FailedToPartitionEngineError', 'GridFireError', 'HashingError', 'IllConditionedJacobianError', 'InvalidQSESolutionError', 'JacobianError', 'KINSolSolverFailureError', 'MissingBaseReactionError', 'MissingKeyReactionError', 'MissingSeedSpeciesError', 'NetworkResizedError', 'PolicyError', 'ReactionError', 'ReactionParsingError', 'SUNDIALSError', 'SingularJacobianError', 'SolverError', 'StaleJacobianError', 'UnableToSetNetworkReactionsError', 'UninitializedJacobianError', 'UnknownJacobianError', 'UtilityError']
class BadCollectionError(EngineError):
    pass
class BadRHSEngineError(EngineError):
    pass
class CVODESolverFailureError(SUNDIALSError):
    pass
class DebugException(GridFireError):
    pass
class EngineError(GridFireError):
    pass
class FailedToPartitionEngineError(EngineError):
    pass
class GridFireError(Exception):
    pass
class HashingError(UtilityError):
    pass
class IllConditionedJacobianError(SolverError):
    pass
class InvalidQSESolutionError(EngineError):
    pass
class JacobianError(EngineError):
    pass
class KINSolSolverFailureError(SUNDIALSError):
    pass
class MissingBaseReactionError(PolicyError):
    pass
class MissingKeyReactionError(PolicyError):
    pass
class MissingSeedSpeciesError(PolicyError):
    pass
class NetworkResizedError(EngineError):
    pass
class PolicyError(GridFireError):
    pass
class ReactionError(GridFireError):
    pass
class ReactionParsingError(ReactionError):
    pass
class SUNDIALSError(SolverError):
    pass
class SingularJacobianError(SolverError):
    pass
class SolverError(GridFireError):
    pass
class StaleJacobianError(JacobianError):
    pass
class UnableToSetNetworkReactionsError(EngineError):
    pass
class UninitializedJacobianError(JacobianError):
    pass
class UnknownJacobianError(JacobianError):
    pass
class UtilityError(GridFireError):
    pass
