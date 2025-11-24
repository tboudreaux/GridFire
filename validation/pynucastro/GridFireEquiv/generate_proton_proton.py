import pynucastro as pyna

from gridfire.engine import GraphEngine, MultiscalePartitioningEngineView, AdaptiveEngineView
from gridfire.type import NetIn

from fourdst.composition import Composition

symbols : list[str] = ["H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"]
X : list[float]     = [0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4]


comp = Composition()
comp.registerSymbol(symbols)
comp.setMassFraction(symbols, X)
comp.finalize(True)

print(f"Initial H-1 mass fraction {comp.getMassFraction("H-1")}")

netIn = NetIn()
netIn.composition = comp
netIn.temperature = 1e7
netIn.density = 1.6e2
netIn.tMax = 1e-9
netIn.dt0 = 1e-12

baseEngine = GraphEngine(netIn.composition, 2)

equiv_species = baseEngine.getNetworkSpecies()
equiv_species = [x.name().replace("-","").lower() for x in equiv_species]

equiv_species = [x if x != 'h1' else 'p' for x in equiv_species]
equiv_species = [x if x != 'n1' else 'n' for x in equiv_species]

rl = pyna.ReacLibLibrary()

# equiv_species = ['p', 'd', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24']
full_lib = rl.linking_nuclei(equiv_species)

print(f"\nFound {len(full_lib.get_rates())} rates.")

net = pyna.PythonNetwork(libraries=[full_lib])

output_filename = "pp_chain_robust.py"
net.write_network(output_filename)

print(f"\nSuccessfully wrote robust PP-chain network to: {output_filename}")
