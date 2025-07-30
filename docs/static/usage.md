# GridFire Python Usage Guide

This tutorial walks you through installing GridFire’s Python bindings, choosing engines and views thoughtfully, running a simulation, and visualizing your results. 

---
## 1. Installation

### 1.1 PyPI Release
The quickest way to get started is:
```bash
pip install gridfire
```

### 1.2 Development from Source
If you want the cutting-edge features or need to hack the C++ backend:
```bash
git clone https://github.com/4DSTAR/GridFire.git
cd GridFire
# Create a virtualenv to isolate dependencies
python3 -m venv .venv && source .venv/bin/activate

# Install Python bindings (meson-python & pybind11 under the hood)
pip install .
```

You can also build manually with Meson (generally end users will not need to do this):
```bash
meson setup build-python
meson compile -C build_gridfire
``` 

---
## 2. Why These Engines and Views?
GridFire’s design [balances physical fidelity and performance](../README.md). Here’s why we pick each component:

1. **GraphEngine**: Constructs the **full** reaction network from Reaclib rates and composition. Use this when:
   - You need maximum physical accuracy (no reactions are omitted).  
   - You are exploring new burning pathways or validating against literature.  

2. **MultiscalePartitioningEngineView**: Implements the Hix & Thielemann partitioning strategy:
   - **Fast reactions** vs **slow reactions** are split onto separate kernels.  
   - This reduces stiffness by isolating processes on very different timescales.  
   - Choose when your network spans orders of magnitude in timescales (e.g., rapid proton captures vs. slow beta decays).  

3. **AdaptiveEngineView**: Dynamically culls low-flow reactions at runtime:
   - At each timestep, reactions with negligible contribution are temporarily removed.  
   - This greatly accelerates large networks without significant loss of accuracy.  
   - Ideal for long integrations where the active set evolves over time.  

4. **Leading-Edge Views**:  
   - `NetworkPrimingEngineView` to inject seed species and study ignition phenomena.  
   - `DefinedEngineView` to freeze the network to a user-specified subset (e.g., focus on the CNO cycle).  

By composing these views in sequence, you can tailor accuracy vs performance for your scientific question. Commonly one might use a flow like **GraphEngine → Partitioning → Adaptive** to capture both full-network physics and manageable stiffness.

---
## 3. Step-by-Step Example
Adapted from [`tests/python/test.py`](../../tests/python/test.py). Comments explain each choice.

```python
import matplotlib.pyplot as plt
from gridfire.engine import GraphEngine, MultiscalePartitioningEngineView, AdaptiveEngineView
from gridfire.solver import DirectNetworkSolver
from gridfire.type import NetIn
from fourdst.composition import Composition

# 1. Define your composition (e.g., M-dwarf surface mix)
symbols = ["H-1","He-3","He-4","C-12","N-14","O-16","Ne-20","Mg-24"]
abundances = [0.708,2.94e-5,0.276,0.003,0.0011,9.62e-3,1.62e-3,5.16e-4]
comp = Composition()
comp.registerSymbols(symbols)
comp.setMassFraction(symbols, abundances)
comp.finalize(normalize=True)  # scale to total mass = 1

# 2. Prepare the NetIn object
netIn = NetIn()
netIn.composition = comp
netIn.temperature = 1.5e7   # Kelvin
netIn.density = 1.6e2       # g/cm³
netIn.tMax = 3.15e7           # seconds to evolve (~1yr)
netIn.dt0 = 1e-12           # initial timestep

# 3. Construct the full network engine
build_depth = 2  # shallow test network for speed; Full = -1
baseEngine = GraphEngine(comp, buildDepth=build_depth)
baseEngine.setUseReverseReactions(False) # At these temps we can turn off reverse reactions

# 4. Partition fast/slow reactions (reduces stiffness)
partitionedEngine = MultiscalePartitioningEngineView(baseEngine)

# 5. Adaptively cull negligible flows (improves speed)
adaptiveEngine = AdaptiveEngineView(partitionedEngine)

# 6. Choose an ODE solver (implicit Rosenbrock4)
solver = DirectNetworkSolver(adaptiveEngine, absTol=1e-12, relTol=1e-8)

# 7. Run the integration

netOut = solver.evaluate(netIn)

# 8. Final result:
print(f"Final H-1 fraction: {netOut.composition.getMassFraction('H-1')}")
```

**Why these choices?**  
- **buildDepth=2**: In Emily’s preliminary tests, depth=2 captures key reaction loops without the overhead of a full network.  
- **Partition & Adaptive Views**: Partitioning reduces stiffness between rapid charged-particle captures and slower β-decays; adaptive culling keeps the working set minimal.  
- **Implicit solver**: Rosenbrock4 handles stiff systems robustly, letting you push to longer `tMax`.

---
## 4. Visualizing Reaction Networks

GridFire engines and views provide built-in export methods for Graphviz DOT and CSV formats:

```python
# Export the base network to DOT for Graphviz
baseEngine.exportToDot('network.dot')
# Optionally generate a PNG (shell):
#   dot -Tpng network.dot -o network.png

# Export a partitioned view of the network
partitionedEngine.exportToDot('partitioned.dot')

# Export to CSV for programmatic analysis
baseEngine.exportToCSV('network.csv')
```  
You can then use tools like Graphviz or pandas:
```bash
# Convert DOT to PNG
dot -Tpng network.dot -o network.png
```  
```python
import pandas as pd

# Load and inspect reaction list
df = pd.read_csv('network.csv')
print(df.head())
```  
For time-series data, record intermediates with an observer and save with pandas or numpy:
```python
import pandas as pd

# Build a DataFrame of time vs species fraction
df = pd.DataFrame({'time': t, 'H-1': X_H1})
df.to_csv('H1_evolution.csv', index=False)
```  
Then plot in pandas or Excel for custom figures.

---
## 5. Beyond the Basics

- **Custom Partition Functions**: In Python, subclass `gridfire.partition.PartitionFunction`, override `evaluate`, `supports`, and `clone` to implement new weighting schemes.  
- **Parameter Studies**: Loop over `buildDepth`, solver tolerances, or initial compositions to get a sense of the sensitity of the network to input conditions or build a monte carlo grid.  
- **Error Handling**:
  ```python
  try:
      results = solver.evaluate(netIn)
  except GridFireRuntimeError as e:
      print('Fatal engine error:', e)
  except GridFireValueError as e:
      print('Invalid input:', e)
  ```

For full API details, consult the docstrings in `src/python/` and the C++ implementation in `src/lib/`. Enjoy exploring nuclear astrophysics with GridFire!
