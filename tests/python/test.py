from gridfire.engine import GraphEngine, MultiscalePartitioningEngineView, AdaptiveEngineView
from gridfire.solver import DirectNetworkSolver
from gridfire.type import NetIn

from fourdst.composition import Composition
from fourdst.atomic import species

symbols : list[str] = ["H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"]
X : list[float]     = [0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4]


comp = Composition()
comp.registerSymbol(symbols)
comp.setMassFraction(symbols, X)
comp.finalize(True)

print(f"Initial H-1 mass fraction {comp.getMassFraction("H-1")}")

netIn = NetIn()
netIn.composition = comp
netIn.temperature = 1.5e7
netIn.density = 1.6e2
netIn.tMax = 4e17
netIn.dt0 = 1e-12

baseEngine = GraphEngine(netIn.composition, 2)
baseEngine.setUseReverseReactions(False)

qseEngine = MultiscalePartitioningEngineView(baseEngine)

adaptiveEngine = AdaptiveEngineView(qseEngine)

solver = DirectNetworkSolver(adaptiveEngine)


def callback(context):
    H1Index = context.engine.getSpeciesIndex(species["H-1"])
    He4Index = context.engine.getSpeciesIndex(species["He-4"])
    C12ndex = context.engine.getSpeciesIndex(species["C-12"])
    Mgh24ndex = context.engine.getSpeciesIndex(species["Mg-24"])
    print(f"Time: {context.t}, H-1: {context.state[H1Index]}, He-4: {context.state[He4Index]}, C-12: {context.state[C12ndex]}, Mg-24: {context.state[Mgh24ndex]}")

# solver.set_callback(callback)
results = solver.evaluate(netIn)

print(f"Final H-1 mass fraction {results.composition.getMassFraction("H-1")}")





