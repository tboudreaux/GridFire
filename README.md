<p align="center">
  <img src="assets/logo/GridFire.png" width="300" alt="OPAT Core Libraries Logo">
</p>


---

GridFire is a C++ library designed to perform general nuclear network evolution using the Reaclib library. It is part of
the larger SERiF project within the 4D-STAR collaboration. GridFire is primarily focused on modeling the most
relevant burning stages for stellar evolution modeling. Currently, there is limited support for inverse reactions.
Therefore, GridFire has a limited set of tools to evolves a fusing plasma in NSE; however,
this is not the primary focus of the library and has therefor not had significant development. For those interested in
modeling super nova, neutron star mergers, or other
high-energy astrophysical phenomena, we **strongly** recomment
using [SkyNet](https://bitbucket.org/jlippuner/skynet/src/master/).

**Design Philosophy and Workflow:**
GridFire is architected to balance physical fidelity, computational efficiency, and extensibility when simulating
complex nuclear reaction networks. Users begin by defining a composition, which is used to construct a full GraphEngine
representation of the reaction network. To manage the inherent stiffness and multiscale nature of these networks,
GridFire employs a layered view strategy: partitioning algorithms isolate fast and slow processes, adaptive culling
removes negligible reactions at runtime, and implicit solvers stably integrate the remaining stiff system. This modular
pipeline allows researchers to tailor accuracy versus performance trade-offs, reuse common engine components, and extend
screening or partitioning models without modifying core integration routines.

## Funding
GridFire is a part of the 4D-STAR collaboration.

4D-STAR is funded by European Research Council (ERC) under the Horizon Europe programme (Synergy Grant agreement No.
101071505: 4D-STAR)
Work for this project is funded by the European Union. Views and opinions expressed are however those of the author(s)
only and do not necessarily reflect those of the European Union or the European Research Council.

## Automatic Build and Installation
### Script Build and Installation Instructions

The easiest way to build GridFire is using the `install.sh` or `install-tui.sh` scripts in the root directory. To use
these scripts, simply run:

```bash
./install.sh
# or
./install-tui.sh
```
The regular installation script will select a standard "ideal" set of build options for you. If you want more control
over the build options, you can use the `install-tui.sh` script, which will provide a text-based user interface to
select the build options you want.

Generally, both are intended to be easy to use and will prompt you automatically to install any missing dependencies.


### Currently known good platforms
The installation script has been tested and found to work on clean installations of the following platforms:
- MacOS 15.3.2 (Apple Silicon + brew installed)
- Fedora 42.0 (aarch64)
- Ubuntu 25.04 (aarch64)
- Ubuntu 22.04 (X86_64)

> **Note:** On Ubuntu 22.04 the user needs to install boost libraries manually as the versions in the Ubuntu repositories
are too old. The installer automatically detects this and will instruct the user in how to do this.

## Manual Build Instructions

### Prerequisites
These only need to be manually installed if the user is not making use of the `install.sh` or `install-tui.sh`

- C++ compiler supporting C++23 standard
- Meson build system (>= 1.5.0)
- Python 3.10 or newer
- CMake 3.20 or newer
- Python packages: `meson-python>=0.15.0`
- Boost libraries (>= 1.75.0) installed system-wide

> **Note:** Boost is the only external library dependency; no additional libraries are required beyond a C++ compiler, Meson, Python, CMake, and Boost.
 
> **Note:** Windows is not supported at this time and *there are no plans to support it in the future*. Windows users are encouraged to use WSL2 or a Linux VM.

### Dependency Installation on Common Platforms

- **Ubuntu/Debian:**
  ```bash
  sudo apt-get update && \
    sudo apt-get install -y build-essential meson python3 python3-pip libboost-all-dev
  ```

- **Fedora/CentOS/RHEL:**
  ```bash
  sudo dnf install -y gcc-c++ meson python3 python3-pip boost-devel
  ```

- **macOS (Homebrew):**
  ```bash
  brew update && \
    brew install boost meson python
  ```

### Building the C++ Library
```bash
meson setup build
meson compile -C build
```

### Installing the Library
```bash
meson install -C build
```

### Python Bindings and Installation
The Python interface is provided via `meson-python` and `pybind11`. To install the Python package:
```bash
pip install .
```

### Developer Workflow
1. Clone the repository and install dependencies listed in `pyproject.toml`.
2. Configure and build with Meson:
   ```bash
   meson setup build
   meson compile -C build
   ```
3. Run the unit tests:
   ```bash
   meson test -C build
   ```
4. Iterate on code, rebuild, and rerun tests.

## Code Architecture and Logical Flow

GridFire is organized into a series of composable modules, each responsible for a specific aspect of nuclear reaction network modeling. The core components include:

- **Engine Module:** Core interfaces and implementations (e.g., `GraphEngine`) that evaluate reaction network rate equations and energy generation.
- **Screening Module:** Implements nuclear reaction screening corrections (`WeakScreening`, `BareScreening`, etc.) affecting reaction rates.
- **Reaction Module:** Parses and manages Reaclib reaction rate data, providing temperature- and density-dependent rate evaluations.
- **Partition Module:** Implements partition functions (e.g., `GroundStatePartitionFunction`, `RauscherThielemannPartitionFunction`) to weight reaction rates based on nuclear properties.
- **Solver Module:** Defines numerical integration strategies (e.g., `DirectNetworkSolver`) for solving the stiff ODE systems arising from reaction networks.
- **Python Interface:** Exposes *almost* all C++ functionality to Python, allowing users to define compositions, configure engines, and run simulations directly from Python scripts.

Generally a user will start by selecting a base engine (currently we only offer `GraphEngine`), which constructs the
full reaction network graph from a given composition. The user can then apply various engine views to adapt the network
topology, such as partitioning fast and slow reactions, adaptively culling low-flow pathways, or priming the network
with specific species. Finally, a numerical solver is selected to integrate the network over time, producing updated
abundances and diagnostics.

### GraphEngine Configuration Options

GraphEngine exposes runtime configuration methods to tailor network construction and rate evaluations:

- **Constructor Parameters:**
  - `BuildDepthType` (`Full`/`Reduced`/`Minimal`): controls network build depth, trading startup time for network completeness.
  - `partition::PartitionFunction`: custom functor for network partitioning based on `Z`, `A`, and `T9`.

- **setPrecomputation(bool precompute):**
  - Enable/disable caching of reaction rates and stoichiometric data at initialization.
  - *Effect:* Reduces per-step overhead; increases memory and setup time.

- **setScreeningModel(ScreeningType type):**
  - Choose plasma screening (models: `BARE`, `WEAK`).
  - *Effect:* Alters rate enhancement under dense/low-T conditions, impacting stiffness.

- **setUseReverseReactions(bool useReverse):**
  - Toggle inclusion of reverse (detailed balance) reactions.
  - *Effect:* Improves equilibrium fidelity; increases network size and stiffness.

### Available Partition Functions

| Function Name                         | Identifier               | Description                                                     |
|---------------------------------------|--------------------------|-----------------------------------------------------------------|
| `GroundStatePartitionFunction`        | "GroundState"          | Weights using nuclear ground-state spin factors.               |
| `RauscherThielemannPartitionFunction` | "RauscherThielemann"   | Interpolates normalized g-factors per Rauscher & Thielemann.   |

These functions implement:
```cpp
double evaluate(int Z, int A, double T9) const;
double evaluateDerivative(int Z, int A, double T9) const;
bool supports(int Z, int A) const;
std::string type() const;
```

## Engine Views

The GridFire engine supports multiple engine view strategies to adapt or restrict network topology. Each view implements a specific algorithm:

| View Name                             | Purpose                                                                          | Algorithm / Reference                                                                                   | When to Use                                                     |
|---------------------------------------|----------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------|
| AdaptiveEngineView                    | Dynamically culls low-flow species and reactions during runtime                 | Iterative flux thresholding to remove reactions below a flow threshold                                   | Large networks to reduce computational cost                     |
| DefinedEngineView                     | Restricts the network to a user-specified subset of species and reactions        | Static network masking based on user-provided species/reaction lists                                      | Targeted pathway studies or code-to-code comparisons            |
| MultiscalePartitioningEngineView      | Partitions the network into fast and slow subsets based on reaction timescales   | Network partitioning following Hix & Thielemann Silicon Burning I & II (DOI:10.1086/177016,10.1086/306692)| Stiff, multi-scale networks requiring tailored integration      |
| NetworkPrimingEngineView              | Primes the network with an initial species or set of species for ignition studies| Single-species injection with transient flow analysis                                                     | Investigations of ignition triggers or initial seed sensitivities|

These engine views implement the common Engine interface and may be composed in any order to build complex network pipelines. New view types can be added by deriving from the `EngineView` base class, and linked into the composition chain without modifying core engine code.

**Python Extensibility:**
Through the Python bindings, users can subclass engine view classes directly in Python, override methods like `evaluate` or `generateStoichiometryMatrix`, and pass instances back into C++ solvers. This enables rapid prototyping of custom view strategies without touching C++ sources.

## Numerical Solver Strategies

GridFire defines a flexible solver architecture through the `networkfire::solver::NetworkSolverStrategy` interface, enabling multiple ODE integration algorithms to be used interchangeably with any engine that implements the `Engine` or `DynamicEngine` contract.

- **NetworkSolverStrategy&lt;EngineT&gt;**: Abstract strategy templated on an engine type. Requires implementation of:
  ```cpp
  NetOut evaluate(const NetIn& netIn);
  ```
  which integrates the network over one timestep and returns updated abundances, temperature, density, and diagnostics.

### DirectNetworkSolver (Implicit Rosenbrock Method)

- **Integrator:** Implicit Rosenbrock4 scheme (order 4) via Boost.Odeint’s `rosenbrock4<double>`, optimized for stiff reaction networks with adaptive step size control using configurable absolute and relative tolerances.
- **Jacobian Assembly:** Employs the `JacobianFunctor` to assemble the Jacobian matrix (∂f/∂Y) at each step, enabling stable implicit integration.
- **RHS Evaluation:** Continues to use the `RHSManager` to compute and cache derivative evaluations and specific energy rates, minimizing redundant computations.
- **Linear Algebra:** Utilizes Boost.uBLAS for state vectors and dense Jacobian matrices, with sparse access patterns supported via coordinate lists of nonzero entries.
- **Error Control and Logging:** Absolute and relative tolerance parameters (`absTol`, `relTol`) are read from configuration; Quill logger captures integration diagnostics and step statistics.

### Algorithmic Workflow in DirectNetworkSolver
1. **Initialization:** Convert input temperature to T9 units, retrieve tolerances, and initialize state vector `Y` from equilibrated composition.
2. **Integrator Setup:** Construct the controlled Rosenbrock4 stepper and bind `RHSManager` and `JacobianFunctor`.
3. **Adaptive Integration Loop:**
   - Perform `integrate_adaptive` advancing until `tMax`, catching any `StaleEngineTrigger` to repartition the network and update composition.
   - On each substep, observe states and log via `RHSManager::observe`.
4. **Finalization:** Assemble final mass fractions, compute accumulated energy, and populate `NetOut` with updated composition and diagnostics.

### Future Solver Implementations
- **Operator Splitting Solvers:** Strategies to decouple thermodynamics, screening, and reaction substeps for performance on stiff, multi-scale networks.
- **GPU-Accelerated Solvers:** Planned use of CUDA/OpenCL backends for large-scale network integration.

These strategies can be developed by inheriting from `NetworkSolverStrategy` and registering against the same engine types without modifying existing engine code.

## Usage Examples

### C++ Example: GraphEngine Initialization
```cpp
#include "gridfire/engine/engine_graph.h"
#include "fourdst/composition/composition.h"

// Define a composition and initialize the engine
fourdst::composition::Composition comp;
gridfire::GraphEngine engine(comp);
```

### C++ Example: Adaptive Network View
```cpp
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/engine/engine_graph.h"

fourdst::composition::Composition comp;
gridfire::GraphEngine baseEngine(comp);
// Dynamically adapt network topology based on reaction flows
gridfire::AdaptiveEngineView adaptiveView(baseEngine);
```

### Python Example
```python
import gridfire
# Initialize GraphEngine with predefined composition
engine = gridfire.GraphEngine(composition="example_composition")
# Perform one integration step
engine.step(dt=1e-3)
print(engine.abundances)
```

## Common Workflow Example

A representative workflow often composes multiple engine views to balance accuracy, stability, and performance when integrating stiff nuclear networks:

```cpp
#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/views/engine_multiscale.h"
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/solver/solver.h"
#include "fourdst/composition/composition.h"

// 1. Define initial composition
fourdst::composition::Composition comp;
// 2. Create base network engine (full reaction graph)
gridfire::GraphEngine baseEngine(comp);

// 3. Partition network into fast/slow subsets (reduces stiffness)
gridfire::MultiscalePartitioningEngineView msView(baseEngine);

// 4. Adaptively cull negligible flux pathways (reduces dimension & stiffness)
gridfire::AdaptiveEngineView adaptView(msView);

// 5. Construct implicit solver (handles remaining stiffness)
gridfire::DirectNetworkSolver solver(adaptView);

// 6. Prepare input conditions
NetIn input{
    comp,     // composition
    1.5e7,      // temperature [K]
    1.5e2,      // density [g/cm^3]
    1e-12,     // initial timestep [s]
    3e17      // integration end time [s]
};

// 7. Execute integration
NetOut output = solver.evaluate(input);
```

**Workflow Components and Effects:**
- **GraphEngine** constructs the full reaction network, capturing all species and reactions.
- **MultiscalePartitioningEngineView** segregates reactions by characteristic timescales (Hix & Thielemann), reducing the effective stiffness by treating fast processes separately.
- **AdaptiveEngineView** prunes low-flux species/reactions at runtime, decreasing dimensionality and improving computational efficiency.
- **DirectNetworkSolver** employs an implicit Rosenbrock method to stably integrate the remaining stiff system with adaptive step control.

This layered approach enhances stability for stiff networks while maintaining accuracy and performance.

## Related Projects

GridFire integrates with and builds upon several key 4D-STAR libraries:

- [fourdst](https://github.com/4D-STAR/fourdst): hub module managing versioning of `libcomposition`, `libconfig`, `liblogging`, and `libconstants`
- [libcomposition](https://github.com/4D-STAR/libcomposition) ([docs](https://4d-star.github.io/libcomposition/)): Composition management toolkit.
- [libconfig](https://github.com/4D-STAR/libconfig): Configuration file parsing utilities.
- [liblogging](https://github.com/4D-STAR/liblogging): Flexible logging framework.
- [libconstants](https://github.com/4D-STAR/libconstants): Physical constants