from gridfire.policy import MainSequencePolicy, NetworkPolicy
from gridfire.engine import DynamicEngine, GraphEngine, EngineTypes
from gridfire.type import NetIn

from typing import Dict
from fourdst.composition import Composition

from testsuite import TestSuite
from utils import init_netIn, init_composition, years_to_seconds

from enum import Enum

EngineNameToType: Dict[str, EngineTypes] = {
    "graphengine": EngineTypes.GRAPH_ENGINE,
    "multiscalepartitioningengineview": EngineTypes.MULTISCALE_PARTITIONING_ENGINE_VIEW,
    "adaptiveengineview": EngineTypes.ADAPTIVE_ENGINE_VIEW
}

class SolarLikeStar_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition()
        super().__init__(
            name="SolarLikeStar_QSE",
            description="GridFire simulation of a roughly solar like star over 10Gyr with QSE enabled.",
            temp=1.5e7,
            density=1.5e2,
            tMax=years_to_seconds(1e10),
            composition=initialComposition,
            notes="Thermodynamically Static, MultiscalePartitioning Engine View"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "AdaptiveEngineView"):
        policy : MainSequencePolicy = MainSequencePolicy(self.composition)
        engine : DynamicEngine = policy.construct()
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()])

class MetalEnhancedSolarLikeStar_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition(ZZs=1)
        super().__init__(
            name="MetalEnhancedSolarLikeStar_QSE",
            description="GridFire simulation of a star with solar core temp and density but enhanced by 1 dex in Z.",
            temp=0.8 * 1.5e7,
            density=1.5e2,
            tMax=years_to_seconds(1e10),
            composition=initialComposition,
            notes="Thermodynamically Static, MultiscalePartitioning Engine View, Z enhanced by 1 dex, temperature reduced to 80% of solar core"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "AdaptiveEngineView"):
        policy : MainSequencePolicy = MainSequencePolicy(self.composition)
        engine : GraphEngine = policy.construct()
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, netIn)
class MetalDepletedSolarLikeStar_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition(ZZs=-1)
        super().__init__(
            name="MetalDepletedSolarLikeStar_QSE",
            description="GridFire simulation of a star with solar core temp and density but depleted by 1 dex in Z.",
            temp=1.2 * 1.5e7,
            density=1.5e2,
            tMax=years_to_seconds(1e10),
            composition=initialComposition,
            notes="Thermodynamically Static, MultiscalePartitioning Engine View, Z depleted by 1 dex, temperature increased to 120% of solar core"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "AdaptiveEngineView"):
        policy : MainSequencePolicy = MainSequencePolicy(self.composition)
        engine : GraphEngine = policy.construct()
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, netIn)

class SolarLikeStar_No_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition()
        super().__init__(
            name="SolarLikeStar_No_QSE",
            description="GridFire simulation of a roughly solar like star over 10Gyr with QSE disabled.",
            temp=1.5e7,
            density=1.5e2,
            tMax=years_to_seconds(1e10),
            composition=initialComposition,
            notes="Thermodynamically Static, No MultiscalePartitioning Engine View"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "AdaptiveEngineView"):
        engine : GraphEngine = GraphEngine(self.composition, 3)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, netIn)


class ValidationSuites(Enum):
    SolarLikeStar_QSE = SolarLikeStar_QSE_Suite
    SolarLikeStar_No_QSE = SolarLikeStar_No_QSE_Suite
    MetalDepletedSolarLikeStar_QSE = MetalDepletedSolarLikeStar_QSE_Suite
    MetalEnhancedSolarLikeStar_QSE = MetalEnhancedSolarLikeStar_QSE_Suite

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run some subset of the GridFire validation suite.")
    parser.add_argument('--suite', type=str, choices=[suite.name for suite in ValidationSuites], nargs="+", help="The validation suite to run.")
    parser.add_argument("--all", action="store_true", help="Run all validation suites.")
    parser.add_argument("--pynucastro-compare", action="store_true", help="Generate pynucastro comparison data.")
    parser.add_argument("--pync-engine", type=str, choices=["GraphEngine", "MultiscalePartitioningEngineView", "AdaptiveEngineView"], default="AdaptiveEngineView", help="The GridFire engine to use to select the reactions for pynucastro comparison.")
    args = parser.parse_args()

    if args.all:
        for suite in ValidationSuites:
            instance : TestSuite = suite.value()
            instance(args.pynucastro_compare, args.pync_engine)

    else:
        for suite_name in args.suite:
            suite = ValidationSuites[suite_name]
            instance : TestSuite = suite.value()
            instance(args.pynucastro_compare, args.pync_engine)

