from gridfire.engine import GraphEngine, MultiscalePartitioningEngineView, AdaptiveEngineView
from gridfire.solver import DirectNetworkSolver
from gridfire.type import NetIn

from fourdst.composition import Composition

symbols : list[str] = ["H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"]
X : list[float]     = [0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4]

comp = Composition();
comp.registerSymbol(symbols)
comp.setMassFraction(symbols, X)
comp.finalize(True)

netIn = NetIn()
netIn.composition = comp
netIn.temperature = 1.5e7
netIn.density = 1.5e2
netIn.tMax = 3e17
netIn.dt0 = 1e-12

baseEngine = GraphEngine(comp, 2)

baseEngine.setUseReverseReactions(False)

qseEngine = MultiscalePartitioningEngineView(baseEngine)
adaptiveEngine = AdaptiveEngineView(qseEngine)

solver = DirectNetworkSolver(adaptiveEngine)
results = solver.evaluate(netIn);
print(results.composition.getMassFraction("H-1"))



