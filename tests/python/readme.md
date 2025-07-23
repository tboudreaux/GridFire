# Python example
GridFire includes robust python bindings. test.py includes a basic example of how to use them.

Generally the concepts are as similar to the C++ API as possible (i.e GridFire code written in C++ and Python will look *very* similar).

```python
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
netIn.tMax = 3e14
netIn.dt0 = 1e-12

baseEngine = GraphEngine(comp, 2)

baseEngine.setUseReverseReactions(False)

qseEngine = MultiscalePartitioningEngineView(baseEngine)
adaptiveEngine = AdaptiveEngineView(qseEngine)

solver = DirectNetworkSolver(adaptiveEngine)
results = solver.evaluate(netIn);
print(results.composition.getMassFraction("H-1"))
```

compare this to  the C++ example and note how similar it is
```c++
#include <iostream>

#include "gridfire/engine/engine.h"
#include "gridfire/solver/solver.h"
#include "gridfire/network.h"
#include "fourdst/composition/composition.h"

int main() {
    using namespace gridfire;
    
    const std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
    const std::vector<double> comp = {0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};


    fourdst::composition::Composition composition;
    composition.registerSymbol(symbols, true);
    composition.setMassFraction(symbols, comp);
    composition.finalize(true);

    NetIn netIn;
    netIn.composition = composition;
    netIn.temperature = 1.5e7;
    netIn.density = 1.5e2;
    netIn.tMax = 3e14;
    netIn.dt0 = 1e-12;

    GraphEngine ReaclibEngine(composition, NetworkBuildDepth::SecondOrder);
    
    ReaclibEngine.setUseReverseReactions(false);
    MultiscalePartitioningEngineView partitioningView(ReaclibEngine);
    AdaptiveEngineView adaptiveView(partitioningView);

    solver::DirectNetworkSolver solver(adaptiveView);
    NetOut netOut;
    netOut = solver.evaluate(netIn);
    std::cout << "H-1 mass fraction: " << netOut.composition.getMassFraction("H-1") << "\n";

}
```

> **Note (1):** The python binding are just that, bindings. They do not include any additional code and simply directly call the underlying C++ code. This means that you inherit almost all of the speed of the C++ code.

## Installing the python bindings
Installing the bindings is as simple as running the following command in the root of the repository:
```bash
pip install .
```
This will use the `mesonpy` build backend to setup and compile the meson project, and then install the python bindings into your current python environment.

> **Note (2):** This will run meson setup and then meson compile to build the project and therefor may take some time and will not produce much output. If you want to see more output during the installation process, you can use the `-v` flag to get verbose output:

If you are a developer and you would like a more convenient edit - compile - test cycle you can install the binding in editable mode
```bash
pip install -e . --no-build-isolation -vv
```
> **Note (3):** Installing in editable mode with `mesonpy` requires that you install the `meson-python' package with pip.

> **Note (4):** The `-vv` flag is optional, but it will give you more verbose output during the installation process, which can be helpful for debugging.

> **Note (5):** The `--no-build-isolation` flag is used to force `mesonpy` to not rebuild the entire project every time and instead to just do incremental builds. meson setup will however **always** run.