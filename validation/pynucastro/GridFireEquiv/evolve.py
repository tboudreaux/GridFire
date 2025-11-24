# import numpy as np
# from scipy.integrate import solve_ivp
#
# import pp_chain_robust as network
#
#
# T = 1e7
# rho = 1.5e2
#
# Y0 = np.zeros(network.nnuc)
# Y0[network.jp] = 0.702583
# Y0[network.jhe3] = 9.74903e-6
# Y0[network.jhe4] = 0.068963
# Y0[network.jc12] = 0.000250029
# Y0[network.jn14] = 7.85632e-5
# Y0[network.jo16] = 0.00060151
# Y0[network.jne20] = 8.10399e-5
# Y0[network.jmg24] = 2.15159e-5
#
# t_start = 0.0
# t_end = 3.14e17
# t_span = (t_start, t_end)
#
# sol = solve_ivp(
#     network.rhs,
#     t_span=t_span,
#     y0=Y0,
#     method='Radau',
#     jac=network.jacobian,
#     args=(rho, T),
#     dense_output=True,
#     rtol=1e-8,
#     atol=1e-8
# )
#
#
# with open("pynucastro_results.csv", 'w') as f:
#     f.write("t,h1,h2,he3,he4,c12,n14,o16,ne20,mg24\n")
#     for (t,h1,h2,he3,he4,c12,n14,o16,ne20,mg24) in zip(sol.t, sol.y[network.jp, :], sol.y[network.jd, :], sol.y[network.jhe3, :], sol.y[network.jhe4, :], sol.y[network.jc12, :], sol.y[network.jn14, :], sol.y[network.jo16, :], sol.y[network.jne20, :], sol.y[network.jmg24, :]):
#         f.write(f"{t},{h1},{h2},{he3},{he4},{c12},{n14},{o16},{ne20},{mg24}\n")
import sys

from gridfire.engine import GraphEngine, MultiscalePartitioningEngineView, AdaptiveEngineView
from gridfire.solver import CVODESolverStrategy
from gridfire.type import NetIn

from fourdst.composition import Composition
from fourdst.atomic import species, Species

from datetime import datetime
from typing import List, Dict, Set, Tuple

X : list[float]     = [0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4]
symbols : list[str] = ["H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"]


comp = Composition()
comp.registerSymbol(symbols)
comp.setMassFraction(symbols, X)
comp.finalize(True)

print(f"Initial H-1 mass fraction {comp.getMassFraction("H-1")}")

netIn = NetIn()
netIn.composition = comp
netIn.temperature = 1.5e7
netIn.density = 1.6e2
netIn.tMax = 3e17
netIn.dt0 = 1e-12

baseEngine = GraphEngine(netIn.composition, 2)
baseEngine.setUseReverseReactions(False)

qseEngine = MultiscalePartitioningEngineView(baseEngine)

adaptiveEngine = AdaptiveEngineView(qseEngine)

solver = CVODESolverStrategy(adaptiveEngine)

data: List[Tuple[float, Dict[str, Tuple[float, float]]]] = []
def callback(context):
    engine = context.engine
    abundances: Dict[str, Tuple[float, float]] = {}
    for species in engine.getNetworkSpecies():
        sid = engine.getSpeciesIndex(species)
        abundances[species.name()] = (species.mass(), context.state[sid])
    data.append((context.t,abundances))

solver.set_callback(callback)

try:
    results = solver.evaluate(netIn, False)
    print(f"H-1 final molar abundance: {results.composition.getMolarAbundance("H-1"):0.3f}")
except:
    print("Warning: solver did not converge", file=sys.stderr)


uniqueSpecies: Set[Tuple[str, float]] = set()
for t, timestep in data:
    for (name, (mass, abundance)) in timestep.items():
        uniqueSpecies.add((name, mass))
sortedUniqueSpecies = list(sorted(uniqueSpecies, key=lambda e: e[1]))

with open(f"gridfire_results_{datetime.now().strftime("%H-%M-%d_%m_%Y")}.csv", 'w') as f:
    f.write('t,')
    for i, (species,mass) in enumerate(sortedUniqueSpecies):
        f.write(f"{species}")
        if i < len(sortedUniqueSpecies)-1:
            f.write(",")
    f.write("\n")

    for t, timestep in data:
        f.write(f"{t},")
        for i, (species, mass) in enumerate(sortedUniqueSpecies):
            if timestep.get(species, None) is not None:
                f.write(f"{timestep.get(species)[1]}")
            else:
                f.write(f"")
            if i < len(sortedUniqueSpecies) - 1:
                f.write(",")
        f.write("\n")

