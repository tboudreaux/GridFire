import os.path
import shutil

from gridfire.policy import MainSequencePolicy, NetworkPolicy
from gridfire.engine import DynamicEngine, GraphEngine, EngineTypes, MultiscalePartitioningEngineView
from gridfire.solver import PointSolverContext
from gridfire.type import NetIn
from gridfire.policy import ConstructionResults

from typing import Dict
from fourdst.composition import Composition

from testsuite import TestSuite
from utils import init_netIn, init_composition, years_to_seconds

from enum import Enum

EngineNameToType: Dict[str, EngineTypes] = {
    "graphengine": EngineTypes.GRAPH_ENGINE,
    "multiscalepartitioningengineview": EngineTypes.MULTISCALE_PARTITIONING_ENGINE_VIEW,
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

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        base_engine: GraphEngine = GraphEngine(self.composition, 3)
        engine : MultiscalePartitioningEngineView = MultiscalePartitioningEngineView(base_engine)
        blob = base_engine.constructStateBlob()
        blob = engine.constructStateBlob(blob)
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

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

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        base_engine: GraphEngine = GraphEngine(self.composition, 3)
        engine : MultiscalePartitioningEngineView = MultiscalePartitioningEngineView(base_engine)
        blob = base_engine.constructStateBlob()
        blob = engine.constructStateBlob(blob)
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

class MetalEnhancedSolarLikeStar_No_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition(ZZs=1)
        super().__init__(
            name="MetalEnhancedSolarLikeStar_No_QSE",
            description="GridFire simulation of a star with solar core temp and density but enhanced by 1 dex in Z.",
            temp=0.8 * 1.5e7,
            density=1.5e2,
            tMax=years_to_seconds(1e10),
            composition=initialComposition,
            notes="Thermodynamically Static, Z enhanced by 1 dex, temperature reduced to 80% of solar core"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        engine: GraphEngine = GraphEngine(self.composition, 4)
        blob = engine.constructStateBlob()
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)


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

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        base_engine: GraphEngine = GraphEngine(self.composition, 3)
        engine : MultiscalePartitioningEngineView = MultiscalePartitioningEngineView(base_engine)
        blob = base_engine.constructStateBlob()
        blob = engine.constructStateBlob(blob)
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

class MetalDepletedSolarLikeStar_No_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition(ZZs=-1)
        super().__init__(
            name="MetalDepletedSolarLikeStar_No_QSE",
            description="GridFire simulation of a star with solar core temp and density but depleted by 1 dex in Z.",
            temp=1.2 * 1.5e7,
            density=1.5e2,
            tMax=years_to_seconds(1e10),
            composition=initialComposition,
            notes="Thermodynamically Static, Z depleted by 1 dex, temperature increased to 120% of solar core"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        engine: GraphEngine = GraphEngine(self.composition, 3)
        blob = engine.constructStateBlob()
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

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

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        engine : GraphEngine = GraphEngine(self.composition, 3)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        context : PointSolverContext = PointSolverContext(engine.constructStateBlob())
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

class SolarLikeStar_No_QSE_Depth_Suite(TestSuite):
    def __init__(self, depth: int = 1):
        initialComposition : Composition = init_composition()
        self.depth : int = depth
        super().__init__(
            name=f"SolarLikeStar_No_QSE_Depth_{depth}_Suite",
            description="GridFire simulation of a roughly solar like star over 10Gyr with QSE disabled.",
            temp=1.5e7,
            density=1.5e2,
            tMax=years_to_seconds(1e10),
            composition=initialComposition,
            notes=f"Thermodynamically Static, No MultiscalePartitioning Engine View, configurable depth {depth}"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        engine : GraphEngine = GraphEngine(self.composition, self.depth)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        context : PointSolverContext = PointSolverContext(engine.constructStateBlob())
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

class HotStar_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition()
        super().__init__(
            name="HotStar_QSE",
            description="GridFire simulation of a hot star over 1Gyr with QSE enabled.",
            temp=2.5e7,
            density=15,
            tMax=1e15,
            composition=initialComposition,
            notes="Thermodynamically Static, MultiscalePartitioning Engine View, B(ish) star conditions"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        base_engine: GraphEngine = GraphEngine(self.composition, 3)
        engine : MultiscalePartitioningEngineView = MultiscalePartitioningEngineView(base_engine)
        blob = base_engine.constructStateBlob()
        blob = engine.constructStateBlob(blob)
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

class CoolStar_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition()
        super().__init__(
            name="CoolStar_QSE",
            description="GridFire simulation of a hot star over 1Gyr with QSE enabled.",
            temp=6e6,
            density=750,
            tMax=1e15,
            composition=initialComposition,
            notes="Thermodynamically Static, MultiscalePartitioning Engine View, M(ish) star conditions"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        base_engine: GraphEngine = GraphEngine(self.composition, 3)
        engine : MultiscalePartitioningEngineView = MultiscalePartitioningEngineView(base_engine)
        blob = base_engine.constructStateBlob()
        blob = engine.constructStateBlob(blob)
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)
class HotStar_No_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition()
        super().__init__(
            name="HotStar_No_QSE",
            description="GridFire simulation of a hot star over 1Gyr with QSE disabled.",
            temp=2.5e7,
            density=15,
            tMax=1e15,
            composition=initialComposition,
            notes="Thermodynamically Static, B(ish) star conditions"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        base_engine: GraphEngine = GraphEngine(self.composition, 3)
        engine : MultiscalePartitioningEngineView = MultiscalePartitioningEngineView(base_engine)
        blob = base_engine.constructStateBlob()
        blob = engine.constructStateBlob(blob)
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

class CoolStar_No_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition()
        super().__init__(
            name="CoolStar_No_QSE",
            description="GridFire simulation of a hot star over 1Gyr with QSE disabled.",
            temp=6e6,
            density=750,
            tMax=1e15,
            composition=initialComposition,
            notes="Thermodynamically Static, M(ish) star conditions"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        base_engine: GraphEngine = GraphEngine(self.composition, 3)
        engine : MultiscalePartitioningEngineView = MultiscalePartitioningEngineView(base_engine)
        blob = base_engine.constructStateBlob()
        blob = engine.constructStateBlob(blob)
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

class VeryCoolStar_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition()
        super().__init__(
            name="VeryCoolStar_QSE",
            description="GridFire simulation of a hot star over 1Gyr with QSE enabled.",
            temp=4e6,
            density=1000,
            tMax=1e15,
            composition=initialComposition,
            notes="Thermodynamically Static, MultiscalePartitioning Engine View, M(ish) star conditions"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        base_engine: GraphEngine = GraphEngine(self.composition, 3)
        engine : MultiscalePartitioningEngineView = MultiscalePartitioningEngineView(base_engine)
        blob = base_engine.constructStateBlob()
        blob = engine.constructStateBlob(blob)
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

class VeryCoolStar_No_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition()
        super().__init__(
            name="VeryCoolStar_No_QSE",
            description="GridFire simulation of a hot star over 1Gyr with QSE disabled.",
            temp=4e6,
            density=1000,
            tMax=1e19,
            composition=initialComposition,
            notes="Thermodynamically Static, M(ish) star conditions"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        engine: GraphEngine = GraphEngine(self.composition, 3)
        blob = engine.constructStateBlob()
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

class VeryHotStar_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition()
        super().__init__(
            name="VeryHotStar_QSE",
            description="GridFire simulation of a hot star over 1Gyr with QSE enabled.",
            temp=4e7,
            density=1,
            tMax=1e15,
            composition=initialComposition,
            notes="Thermodynamically Static, MultiscalePartitioning Engine View, M(ish) star conditions"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        base_engine: GraphEngine = GraphEngine(self.composition, 4)
        engine : MultiscalePartitioningEngineView = MultiscalePartitioningEngineView(base_engine)
        blob = base_engine.constructStateBlob()
        blob = engine.constructStateBlob(blob)
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)

class VeryHotStar_No_QSE_Suite(TestSuite):
    def __init__(self):
        initialComposition : Composition = init_composition()
        super().__init__(
            name="VeryHotStar_No_QSE",
            description="GridFire simulation of a hot star over 1Gyr with QSE disabled.",
            temp=4e7,
            density=1,
            tMax=1e15,
            composition=initialComposition,
            notes="Thermodynamically Static, B(ish) star conditions"
        )

    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        engine: GraphEngine = GraphEngine(self.composition, 4)
        blob = engine.constructStateBlob()
        context : PointSolverContext = PointSolverContext(blob)
        netIn : NetIn = init_netIn(self.temperature, self.density, self.tMax, self.composition)
        self.evolve(engine, context, netIn, pynucastro_compare = pynucastro_compare, engine_type=EngineNameToType[pync_engine.lower()], output=output)


class ValidationSuites(Enum):
    SolarLikeStar_QSE = SolarLikeStar_QSE_Suite
    SolarLikeStar_No_QSE = SolarLikeStar_No_QSE_Suite
    SolarLikeStar_No_QSE_Depth = SolarLikeStar_No_QSE_Depth_Suite
    MetalDepletedSolarLikeStar_QSE = MetalDepletedSolarLikeStar_QSE_Suite
    MetalEnhancedSolarLikeStar_QSE = MetalEnhancedSolarLikeStar_QSE_Suite
    MetalDepletedSolarLikeStar_No_QSE = MetalDepletedSolarLikeStar_No_QSE_Suite
    MetalEnhancedSolarLikeStar_No_QSE = MetalEnhancedSolarLikeStar_No_QSE_Suite
    HotStar_QSE = HotStar_QSE_Suite
    CoolStar_QSE = CoolStar_QSE_Suite
    HotStar_No_QSE = HotStar_No_QSE_Suite
    CoolStar_No_QSE = CoolStar_No_QSE_Suite
    VeryCoolStar_QSE = VeryCoolStar_QSE_Suite
    VeryHotStar_QSE = VeryHotStar_QSE_Suite
    VeryCoolStar_No_QSE = VeryCoolStar_No_QSE_Suite
    VeryHotStar_No_QSE = VeryHotStar_No_QSE_Suite

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run some subset of the GridFire validation suite.")
    parser.add_argument('--suite', type=str, choices=[suite.name for suite in ValidationSuites], nargs="+", help="The validation suite to run.")
    parser.add_argument("--all", action="store_true", help="Run all validation suites.")
    parser.add_argument("--pynucastro-compare", action="store_true", help="Generate pynucastro comparison data.")
    parser.add_argument("--pync-engine", type=str, choices=["GraphEngine", "MultiscalePartitioningEngineView", "AdaptiveEngineView"], default="AdaptiveEngineView", help="The GridFire engine to use to select the reactions for pynucastro comparison.")
    parser.add_argument("-o", "--output", type=str, help="Directory to save OKAY results too", default="GF_Validation_Output")
    parser.add_argument("--depths", type=int, nargs="+", default=[1, 2, 3, 4, 5, 6, 7], help="Construction depths to test. Must be positive non zero")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    if args.all:
        for suite in ValidationSuites:
            if suite.name == "SolarLikeStar_No_QSE_Depth":
                    for depth in args.depths:
                        instance : TestSuite = suite.value(depth=depth)
                        instance(args.pynucastro_compare, args.pync_engine, args.output)
            else:
                instance : TestSuite = suite.value()
                instance(args.pynucastro_compare, args.pync_engine, args.output)

    else:
        for suite_name in args.suite:
            suite = ValidationSuites[suite_name]
            if suite.name == "SolarLikeStar_No_QSE_Depth":
                    for depth in args.depths:
                        instance : TestSuite = suite.value(depth=depth)
                        instance(args.pynucastro_compare, args.pync_engine, args.output)
            else:
                instance : TestSuite = suite.value()
                instance(args.pynucastro_compare, args.pync_engine, args.output)

