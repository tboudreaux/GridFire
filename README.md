<p align="center">
  <img src="https://github.com/4D-STAR/GridFire/blob/main/assets/logo/GridFire.png?raw=true" width="300" alt="OPAT Core Libraries Logo">
</p>

---
![PyPI - Version](https://img.shields.io/pypi/v/gridfire?style=for-the-badge)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/gridfire?style=for-the-badge)

![GitHub License](https://img.shields.io/github/license/4D-STAR/GridFire?style=for-the-badge)
![ERC](https://img.shields.io/badge/Funded%20by-ERC-blue?style=for-the-badge&logo=europeancommission)

![Dynamic Regex Badge](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fgithub.com%2F4D-STAR%2FGridFire%2Fblob%2Fmain%2Fmeson.build&search=version%3A%20'(%5B0-9a-zA-Z%5C.%5D%2B)'&style=for-the-badge&label=GitHub%20Main%20Branch)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/4D-STAR/GridFire?style=for-the-badge)


---

# Table of Contents
- [Introduction](#introduction)
  - [Design Philosophy and Workflow](#design-philosophy-and-workflow)
  - [Funding](#funding)
- [Usage](#usage)
  - [Python installation](#python-installation)
    - [pypi](#pypi)
    - [source](#source)
    - [source for developers](#source-for-developers)
    - [patching shared object files](#patching-shared-object-files)
  - [Automatic Build and Installation](#automatic-build-and-installation)
    - [Script Build and Installation Instructions](#script-build-and-installation-instructions)
    - [Currently, known good platforms](#currently-known-good-platforms)
  - [Manual Build Instructions](#manual-build-instructions)
    - [Prerequisites](#prerequisites)
    - [Install Scripts](#install-scripts)
    - [Dependency Installation on Common Platforms](#dependency-installation-on-common-platforms)
    - [Building the C++ Library](#building-the-c-library)
    - [Installing the Library](#installing-the-library)
    - [Minimum compiler versions](#minimum-compiler-versions)
- [Code Architecture and Logical Flow](#code-architecture-and-logical-flow)
- [Engines](#engines)
  - [GraphEngine](#graphengine)
  - [GraphEngine Configuration Options](#graphengine-configuration-options)
  - [Available Partition Functions](#available-partition-functions)
  - [AutoDiff](#autodiff)
- [Reaclib in GridFire](#reaclib-in-gridfire)
- [Engine Views](#engine-views)
  - [A Note about composability](#a-note-about-composability)
- [Numerical Solver Strategies](#numerical-solver-strategies)
  - [NetworkSolverStrategy&lt;EngineT&gt;](#networksolverstrategyenginet)
  - [NetIn and NetOut](#netin-and-netout)
  - [DirectNetworkSolver (Implicit Rosenbrock Method)](#directnetworksolver-implicit-rosenbrock-method)
  - [Algorithmic Workflow in DirectNetworkSolver](#algorithmic-workflow-in-directnetworksolver)
  - [Future Solver Implementations](#future-solver-implementations)
- [Python Extensibility](#python-extensibility)
- [Usage Examples](#usage-examples)
  - [C++](#c)
    - [GraphEngine Initialization](#graphengine-initialization)
    - [Adaptive Network View](#adaptive-network-view)
    - [Composition Initialization](#composition-initialization)
    - [Common Workflow Example](#common-workflow-example)
    - [Callback Example](#callback-example)
  - [Python](#python)
    - [Common Workflow Example](#common-workflow-example-1)
    - [Python callbacks](#python-callbacks)
- [Related Projects](#related-projects)

# Introduction
GridFire is a C++ library designed to perform general nuclear network
evolution. It is part of the larger SERiF project within the 4D-STAR
collaboration. GridFire is primarily focused on modeling the most relevant
burning stages for stellar evolution modeling. Currently, there is limited
support for inverse reactions. Therefore, GridFire has a limited set of tools
to evolve a fusing plasma in NSE; however, this is not the primary focus of
the library and has therefore not had significant development. For those
interested in modeling super nova, neutron star mergers, or other high-energy
astrophysical phenomena, we **strongly** recommend using
[SkyNet](https://bitbucket.org/jlippuner/skynet/src/master/).

## Design Philosophy and Workflow
GridFire is architected to balance physical fidelity, computational efficiency,
and extensibility when simulating complex nuclear reaction networks. Users
begin by defining a composition, which is used to construct a full GraphEngine
representation of the reaction network. A GraphNetwork uses [JINA
Reaclib](https://reaclib.jinaweb.org/index.php) reaction rates ([Cyburt et al., ApJS 189
(2010) 240.](https://iopscience.iop.org/article/10.1088/0067-0049/189/1/240)) along
with a dynamically constructed network topology. To manage the inherent
stiffness and multiscale nature of these networks, GridFire employs a layered
view strategy: partitioning algorithms isolate fast and slow processes,
adaptive culling removes negligible reactions at runtime, and implicit solvers
stably integrate the remaining stiff system. This modular pipeline allows
researchers to tailor accuracy versus performance trade-offs, reuse common
engine components, and extend screening or partitioning models without
modifying core integration routines.

## Funding
GridFire is a part of the 4D-STAR collaboration.

4D-STAR is funded by European Research Council (ERC) under the Horizon Europe
programme (Synergy Grant agreement No. 101071505: 4D-STAR) Work for this
project is funded by the European Union. Views and opinions expressed are
however those of the author(s) only and do not necessarily reflect those of the
European Union or the European Research Council.

# Usage

## Python installation
By far the easiest way to install is with pip. This will install either
pre-compiled wheels or, if your system has not had a wheel compiled for it, it
will try to build locally (this may take **a long time**). The python bindings
are just that and should maintain nearly the same speed as the C++ code. End
users are strongly encouraged to use the python module rather than the C++ code.

### pypi
Installing from pip is as simple as.
```bash
pip install gridfire
```

These wheels have been compiled on many systems

| Version   | Platform | Architecture | CPython Versions                                           | PyPy Versions |
|-----------|----------|--------------|------------------------------------------------------------|---------------|
 | 0.7.0_rc1 | macOS    | arm64        | 3.8, 3.9, 3.10, 3.11, 3.12, 3.13 (std & t), 3.14 (std & t) | 3.10, 3.11    |
| 0.7.0_rc1 | Linux    | aarch64      | 3.8, 3.9, 3.10, 3.11, 3.12, 3.13 (std & t), 3.14 (std & t) | 3.10, 3.11    |
| 0.7.0_rc1 | Linux    | x86\_64      | 3.8, 3.9, 3.10, 3.11, 3.12, 3.13 (std & t), 3.14 (std & t) | 3.10, 3.11    |
| 0.5.0     | macOS    | arm64        | 3.8, 3.9, 3.10, 3.11, 3.12, 3.13 (std & t), 3.14 (std & t) | 3.10, 3.11    |
| 0.5.0     | Linux    | aarch64      | 3.8, 3.9, 3.10, 3.11, 3.12, 3.13 (std & t), 3.14 (std & t) | 3.10, 3.11    |
| 0.5.0     | Linux    | x86\_64      | 3.8, 3.9, 3.10, 3.11, 3.12, 3.13 (std & t), 3.14 (std & t) | 3.10, 3.11    |

> **Note**: Currently macOS x86\_64 does **not** have a precompiled wheel. Due
> to that platform being phased out it is likely that there will never be
> precompiled wheels or releases for it.

> **Note:** macOS wheels were targeted to macOS 12 Monterey and should work on
> any version more recent than that (at least as of August 2025).

> **Note:** Linux wheels were compiled using manylinux_2_28 and are expected to
> work on Debian 10+, Ubuntu 18.10+, Fedora 29+, or CentOS/RHEL 8+

> **Note:** If your system does not have a prebuilt wheel the source
> distribution will download from pypi and try to build. This may simply not
> work if you do not have the correct system dependencies installed. If it
> fails the best bet is to try to build boost >= 1.83.0 from source and install
> (https://www.boost.org/) as that is the most common broken dependency.

### source
The user may also build the python bindings directly from source

```bash
git clone https://github.com/4D-STAR/GridFire
cd GridFire
pip install .
```

> **Note:** that if you do not have all system dependencies installed this will
> fail, the steps in further sections address these in more detail.

> **Note:** If you are using macos you should use the included `pip_install_mac_patch.sh` script instead of `pip install .` as this will automatically patch the build shared object libraries such that they can be loaded by the macos dynamic loader.

### source for developers
If you are a developer and would like an editable and incremental python
install `meson-python` makes this very easy

```bash
git clone https://github.com/4D-STAR/GridFire
cd GridFire
pip install -e . --no-build-isolation -vv
```

This will generate incremental builds whenever source code changes, and you run
a python script automatically (note that since `meson setup` must run for each
of these it does still take a few seconds to recompile regardless of how small
a source code change you have made). It is **strongly** recommended that
developers use this approach and end users *do not*.

#### Patching Shared Object Files
If you need to patch shared object files generated by meson-python directly you should first locate the shared object file
these will be in the site-packages and site-packages/fourdst directories for your python environment.

Look for files named

- `site-packages/gridfire.cpython-3*-darwin.so`
- `site-packages/fourdst/_phys.cpython-3*-darwin.so` 

then, for each of these files, run

```bash
otool -l <Path/to/file> | grep RPATH -A2
```

count the number of occurrences of duplicate RPATH entries (these should look like `@loaderpath/.gridfire.mesonpy.libs` or `@loaderpath/../.fourdst.mesonpy.libs`). Then use `install_name_tool` to remove **all but one of these** from each shared object file.

If for example there are 4 occurrences of the path `@loader_path/../.fourdst.mesonpy.libs` in `_phys.cpython-3*-darwin.so` then you should run the following command 3 times
```bash
install_name_tool -delete_rpath @loader_path/../.fourdst.mesonpy.libs site-packages/fourdst/_phys.cpython-314-darwin.so
```

the same for the other shared object file (make sure to count the duplicate rpath entries for each separately as there may be a different number of duplicates in each shared object file).

We also include a script at `pip_install_mac_patch.sh` which will do this automatically for you.

## Automatic Build and Installation
### Script Build and Installation Instructions

The easiest way to build GridFire is using the `install.sh` or `install-tui.sh`
scripts in the root directory. To use these scripts, simply run:

```bash
./install.sh
# or
./install-tui.sh
```
The regular installation script will select a standard "ideal" set of build
options for you. If you want more control over the build options, you can use
the `install-tui.sh` script, which will provide a text-based user interface to
select the build options you want.

Generally, both are intended to be easy to use and will prompt you
automatically to install any missing dependencies.


### Currently, known good platforms
The installation script has been tested and found to work on clean
installations of the following platforms:
- MacOS 15.3.2 (Apple Silicon + brew installed)
- Fedora 42.0 (aarch64)
- Ubuntu 25.04 (aarch64)
- Ubuntu 22.04 (X86_64)

> **Note:** On Ubuntu 22.04 the user needs to install boost libraries manually
> as the versions in the Ubuntu repositories
> are too old. The installer automatically detects this and will instruct the
> user in how to do this.

## Manual Build Instructions

### Prerequisites
These only need to be manually installed if the user is not making use of the
`install.sh` or `install-tui.sh`

#### Required
- C++ compiler supporting C++23 standard
- Meson build system (>= 1.5.0)
- Python 3.8 or newer
- CMake 3.20 or newer
- ninja 1.10.0 or newer
- Python packages: `meson-python>=0.15.0`
- Boost libraries (>= 1.83.0) installed system-wide (or at least findable by
  meson with pkg-config)

#### Optional
- dialog (used by the `install.sh` script, not needed if using pip or meson
  directly)
- pip (used by the `install.sh` script or by calling pip directly, not needed
  if using meson directly)

> **Note:** Boost is the only external library dependency used by GridFire directly.

> **Note:** Windows is not supported at this time and *there are no plans to
> support it in the future*. Windows users are encouraged to use WSL2 or a
> Linux VM.

> **Note:** If `install-tui.sh` is not able to find a usable version of boost
> it will provide directions to fetch, compile, and install a usable version.

### Install Scripts
GridFire ships with an installer (`install.sh`) which is intended to make the
process of installation both easier and more repeatable.

#### Ease of Installation
Both scripts are intended to automate installation more or less completely. This
includes dependency checking. In the event that a dependency cannot be found
they try to install (after explicitly asking for user permission). If that does
not work they will provide a clear message as to what went wrong.

#### Reproducibility
The TUI mode provides easy modification of meson build system and compiler
settings which can then be saved to a config file. This config file can then be
loaded by either tui mode or cli mode (with the `--config`) flag meaning that
build configurations can be made and reused. Note that this is **not** a
deterministically reproducible build system as it does not interact with any
system dependencies or settings, only meson and compiler settings.

#### Examples

##### TUI config and saving
[![asciicast](https://asciinema.org/a/ahIrQPL71ErZv5EKKujfO1ZEW.svg)](https://asciinema.org/a/ahIrQPL71ErZv5EKKujfO1ZEW)

##### TUI config loading and meson setup
[![asciicast](https://asciinema.org/a/zGdzt9kYsETltG0TJKC50g3BK.svg)](https://asciinema.org/a/zGdzt9kYsETltG0TJKC50g3BK)

##### CLI config loading, setup, and build
[![asciicast](https://asciinema.org/a/GYaWTXZbDJRD4ohde0s3DkFMC.svg)](https://asciinema.org/a/GYaWTXZbDJRD4ohde0s3DkFMC)


> **Note:** `install-tui.sh` is simply a script which calls `install.sh` with
> the `--tui` flag. You can get the exact same results by running `install.sh
> --tui`.

> **Note:** Call `install.sh` with the `--help` or `--h` flag to see command
> line options

> **Note:** `clang` tends to compile GridFire much faster than `gcc` thus why I
> select it in the above asciinema recording.

### Dependency Installation on Common Platforms

- **Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install -y build-essential meson python3 python3-pip libboost-all-dev
```

> **Note:** Depending on the ubuntu version you have the libboost-all-dev
> libraries may be too old. If this is the case refer to the boost
> documentation for how to download and install a version `>=1.83.0`

> **Note:** On recent versions of ubuntu python has switched to being
> externally managed by the system. We **strongly** recommend that if you
> install manually all python packages are installed inside some kind of
> virtual environment (e.g. `pyenv`, `conda`, `python-venv`, etc...). When using
> the installer script this is handled automatically using `python-venv`.

- **Fedora/CentOS/RHEL:**
```bash
sudo dnf install -y gcc-c++ meson python3 python3-pip boost-devel
```

- **macOS (Homebrew):**
```bash
brew update
brew install boost meson python
```

### Building the C++ Library
```bash
meson setup build
meson compile -C build
```

#### Clang vs. GCC
As noted above `clang` tends to compile GridFire much faster than `gcc`. If
your system has both `clang` and `gcc` installed you may force meson to use
clang via environmental variables

```bash
CC=clang CXX=clang++ meson setup build_clang
meson compile -C build_clang
```

### Installing the Library
```bash
meson install -C build
```

### Minimum compiler versions
GridFire uses C++23 features and therefore only compilers and standard library
implementations which support C++23 are supported. Generally we have found that
`gcc >= 13.0.0` or `clang >= 16.0.0` work well.


## Code Architecture and Logical Flow

GridFire is organized into a series of composable modules, each responsible for
a specific aspect of nuclear reaction network modeling. The core components
include:

- **Engine Module:** Core interfaces and implementations (e.g., `GraphEngine`)
  that evaluate reaction network rate equations and energy generation. Also
  implemented `Views` submodule.
- **Engine::Views Module:** Composable engine optimization and modification
  (e.g. `MultiscalePartitioningEngineView`) which can be used to make a problem
  more tractable or applicable.
- **Screening Module:** Implements nuclear reaction screening corrections (e.g.
  `WeakScreening` ([Salpeter,
  1954](https://adsabs.harvard.edu/full/1954AuJPh...7..373S)), `BareScreening`)
  affecting reaction rates.
- **Reaction Module:** Parses and manages Reaclib reaction rate data, providing
  temperature- and density-dependent rate evaluations.
- **Partition Module:** Implements partition functions (e.g.,
  `GroundStatePartitionFunction`, `RauscherThielemannPartitionFunction`
  ([Rauscher & Thielemann,
  2000](https://www.sciencedirect.com/science/article/pii/S0092640X00908349?via%3Dihub]))
  to weight reaction rates based on nuclear properties.
- **Solver Module:** Defines numerical integration strategies (e.g.,
  `CVODESolverStrategy`) for solving the stiff ODE systems arising from reaction
  networks.
- **io Module:** Defines shared interface for parsing network data from files
- **trigger Module:** Defines interface for complex trigger logic so that repartitioning can be followed.
- **Policy Module:** Contains "policies" which are small modular units of code that enforce certain contracts. 
For example the `ProtonProtonReactionChainPolicy` enforces than an engine must include at least all the reactions
in the proton-proton chain. This module exposes the primary construction interface for users. I.e. select a policy (such as `MainSequencePolicy`), provide a composition, and get back an engine which satisfies that policy.
- **Python Interface:** Exposes *almost* all C++ functionality to Python,
  allowing users to define compositions, configure engines, and run simulations
  directly from Python scripts.

Generally a user will start by selecting a base engine (currently we only offer
`GraphEngine`), which constructs the full reaction network graph from a given
composition. The user can then apply various engine views to adapt the network
topology, such as partitioning fast and slow reactions, adaptively culling
low-flow pathways, or priming the network with specific species. Finally, a
numerical solver is selected to integrate the network over time, producing
updated
abundances and diagnostics.

## Engines
GridFire is, at its core, based on a series of `Engines`. These are constructs
which know how to report information on series of ODEs which need to be solved
to evolve abundances. The important thing to understand about `Engines` is
that they contain all the detailed physics GridFire uses. For example a
`Solver` takes an `Engine` but does not compute physics itself. Rather, it asks
the `Engine` for stuff like the jacobian matrix, stoichiometry, nuclear energy
generation rate, and change in abundance with time.

Refer to the API documentation for the exact interface which an `Engine` must
implement to be compatible with GridFire solvers.

Currently, we only implement `GraphEngine` which is intended to be a very general and
adaptable `Engine`.

### GraphEngine
In GridFire the `GraphEngine` will generally be the most fundamental building
block of a nuclear network. A `GraphEngine` represents a directional hypergraph
connecting some set of atomic species through reactions listed in the [JINA
Reaclib database](https://reaclib.jinaweb.org/index.php).

`GraphEngine`s are constructed from a seed composition of species from which
they recursively expand their topology outward, following known reaction
pathways and adding new species to the tracked list as they expand.


### GraphEngine Configuration Options

GraphEngine exposes runtime configuration methods to tailor network
construction and rate evaluations:

- **Constructor Parameters:**
    - `composition`: The initial seed composition to start network construction from.
    - `BuildDepthType` (`Full`, `Shallow`, `SecondOrder`, etc...): controls
      number of recursions used to construct the network topology. Can either be a
      member of the `NetworkBuildDepth` enum or an integer.
    - `partition::PartitionFunction`: Partition function used when evaluating
      detailed balance for inverse rates.
    - `NetworkConstructionFlags`: A bitwise flag telling the network how to construct itself. That is, what reaction types should be used in construction. For example one might use `NetworkConstructionFlags::STRONG | NetworkConstructionFlags::BETA_PLUS` to use all strong reactions and β+ decay. By Default this is set to use reaclib strong and reaclib weak (no WRL included by default due to current pathological stiffness issues).

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

| Function Name                         | Identifier / Enum    | Description                                                                                                              |
|---------------------------------------|----------------------|--------------------------------------------------------------------------------------------------------------------------|
| `GroundStatePartitionFunction`        | "GroundState"        | Weights using nuclear ground-state spin factors.                                                                         |
| `RauscherThielemannPartitionFunction` | "RauscherThielemann" | Interpolates normalized g-factors per Rauscher & Thielemann.                                                             |
| `CompositePartitionFunction`          | "Composite"          | Combines multiple partition functions for situations where different partitions functions are used for different domains |

### AutoDiff
One of the primary tasks any engine must accomplish is to report the jacobian
matrix of the system to the solver. `GraphEngine`  uses `CppAD`, a C++ auto
differentiation library, to generate analytic jacobian matrices very
efficiently.


## Reaclib in GridFire
All reactions in JINA Reaclib which only include reactants iron and lighter
were downloaded on June 17th, 2025 where the most recent documented change on
the JINA Reaclib site was on June 24th, 2021.

All of these reactions have been compiled into a header file which is then
statically compiled into the gridfire binaries (specifically into
lib_reaction_reaclib.cpp.o). This does increase the binary size by a few MB;
however, the benefit is faster load times and more importantly no need for end
users to manage resource files.

If a developer wants to add new reaclib reactions we include a script at
`utils/reaclib/format.py` which can ingest a reaclib data file and produce the
needed header file. More details on this process are included in
`utils/reaclib/readme.md`


## Engine Views

The GridFire engine supports multiple engine view strategies to adapt or
restrict network topology. Generally when extending GridFire the approach is
likely to be one of adding new `EngineViews`.

| View Name                        | Purpose                                                                           | Algorithm / Reference                                                                                      | When to Use                                                       |
|----------------------------------|-----------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------|
| AdaptiveEngineView               | Dynamically culls low-flow species and reactions during runtime                   | Iterative flux thresholding to remove reactions below a flow threshold                                     | Large networks to reduce computational cost                       |
| DefinedEngineView                | Restricts the network to a user-specified subset of species and reactions         | Static network masking based on user-provided species/reaction lists                                       | Targeted pathway studies or code-to-code comparisons              |
| FileDefinedEngineView            | Load a defined engine view from a file using some parser                          | Same as DefinedEngineView but loads from a file                                                            | Same as DefinedEngineView                                         |
| MultiscalePartitioningEngineView | Partitions the network into fast and slow subsets based on reaction timescales    | Network partitioning following Hix & Thielemann Silicon Burning I & II (DOI:10.1086/177016,10.1086/306692) | Stiff, multi-scale networks requiring tailored integration        |
| NetworkPrimingEngineView         | Primes the network with an initial species or set of species for ignition studies | Single-species ignition and network priming                                                                | Investigations of ignition triggers or initial seed sensitivities |

These engine views implement the common Engine interface and may be composed in
any order to build complex network pipelines. New view types can be added by
deriving from the `EngineView` base class, and linked into the composition
chain without modifying core engine code.

### A Note about composability
There are certain functions for which it is expected that a call to an engine
view will propagate the result down the chain of engine views, eventually
reaching the base engine (e.g. `DynamicEngine::update`). We do not strongly
enforce this as it is not hard to contrive a situation where that is not the
mose useful behavior; however, we do strongly encourage developers to think
carefully about passing along calls to base engine methods when implementing
new views.

## Numerical Solver Strategies

GridFire defines a flexible solver architecture through the
`networkfire::solver::NetworkSolverStrategy` interface, enabling multiple ODE
integration algorithms to be used interchangeably with any engine that
implements the `Engine` or `DynamicEngine` contract.

### NetworkSolverStrategy&lt;EngineT&gt;:
All GridFire solvers implement the abstract strategy templated by
`NetworkSolverStrategy` which enforces only that there is some `evaluate`
method with the following signature

```cpp
NetOut evaluate(const NetIn& netIn);
```
Which is intended to integrate some network over some time and returns updated
abundances, temperature, density, and diagnostics.

### NetIn and NetOut
GridFire solvers use a unified input and output type for their public interface
(though as developers will quickly learn, internally these are immediately
broken down into simpler data structures). All solvers expect a `NetIn` struct
for the input type to the `evaluate` method and return a `NetOut` struct.

#### NetIn
A `NetIn` struct contains
- The composition to start the timestep at. (`NetIn::composition`)
- The temperature in Kelvin (`NetIn::temperature`)
- The density in g/cm^3 (`NetIn::density`)
- The max time to evolve the network to in seconds (`NetIn::tMax`)
- The initial timestep to use in seconds (`NetIn::dt0`)
- The initial energy in the system in ergs (`NetIn::energy`)

>**Note:** It is often useful to set `NetIn::dt0` to something *very* small and
>let an iterative time stepper push the timestep up. Often for main sequence
>burning I use ~1e-12 for dt0

>**Note:** The composition must be a `fourdst::composition::Composition`
>object. This is made available through the `foursdt` library and the
>`fourdst/composition/Composition.h` header. `fourdst` is installed
>automatically with GridFire

>**Note:** In Python composition comes from `fourdst.composition.Composition`
>and similarly is installed automatically when building GridFire python
>bindings.

#### NetOut
A `NetOut` struct contains
- The final composition after evolving to `tMax` (`NetOut::composition`)
- The number of steps the solver took to evolve to `tmax` (`NetOut::num_steps`)
- The final energy generated by the network while evolving to `tMax`
  (`NetOut::energy`)

>**Note:** Currently `GraphEngine` only considers energy due to nuclear mass
>defect and not neutrino loss.

### CVODESolverStrategy

We use the CVODE module from [SUNDIALS](https://computing.llnl.gov/projects/sundials/cvode) as our primary numerical
solver. Specifically we use the BDF linear multistep method from that which includes advanced adaptive timestepping.

Further, we use a trigger system to periodically repartition the network as the state of the network changes. This
keeps the stiffness of the network tractable. The algorithm we use for that is 

1. Trigger every 1000th time that the simulation time exceeds the simulationTimeInterval
2. OR if any off-diagonal Jacobian entry exceeds the offDiagonalThreshold
3. OR every 10th time that the timestep growth exceeds the timestepGrowthThreshold (relative or absolute)
4. OR if the number of convergence failures grows more than 100% from one step to the next or exceeds 5 at any given step.

Moreover, callback functions can be registered in either python or C++ which will take a `const CVODESolverStrategy::TimestepContext&` struct 
as argument. This allows for more complex logging logic. Note that callbacks **do not** let you reach inside the 
solver and adjust the state of the network. They are only intended for investigation not extension of physics. If you
wish to extend the physics this must be implemented at the engine or engine view level.

## Python Extensibility
Through the Python bindings, users can subclass engine view classes directly in
Python, override methods like `evaluate` or `generateStoichiometryMatrix`, and
pass instances back into C++ solvers. This enables rapid prototyping of custom
view strategies without touching C++ sources.


# Usage Examples

## C++

### GraphEngine Initialization
```c++
#include "gridfire/engine/engine_graph.h"
#include "fourdst/composition/composition.h"

int main(){
    // Define a composition and initialize the engine
    fourdst::composition::Composition comp;
    gridfire::GraphEngine engine(comp);
}
```

### Adaptive Network View
```c++
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/engine/engine_graph.h"

int main(){
    fourdst::composition::Composition comp;
    gridfire::GraphEngine baseEngine(comp);
    // Dynamically adapt network topology based on reaction flows
    gridfire::AdaptiveEngineView adaptiveView(baseEngine);
}
```

### Composition Initialization
```c++
#include "fourdst/composition/composition.h"
#include "fourdst/composition/utils.h" // for buildCompositionFromMassFractions
#include <vector>
#include <string>

#include <iostream>

int main() {

    std::vector<std::string> symbols = {"H-1", "He-4", "C-12"};
    std::vector<double> massFractions = {0.7, 0.29, 0.01};
    
    const fourdst::composition::Composition comp = fourdst::composition::buildCompositionFromMassFractions(symbols, massFractions);

    std::cout << comp << std::endl;
}
```

### Common Workflow Example

A representative workflow often composes multiple engine views to balance
accuracy, stability, and performance when integrating stiff nuclear networks:

```c++
#include "gridfire/gridfire.h" // Unified header for real usage

#include "fourdst/composition/composition.h"
#include "fourdst/composition/utils.h" // for buildCompositionFromMassFractions

int main(){
    // 1. Define initial composition
    std::unordered_map<std::string, double> initialMassFractions = {
        {"H-1", 0.7},
        {"He-4", 0.29},
        {"C-12", 0.01}
    };
    const fourdst::composition::Composition composition = fourdst::composition::buildCompositionFromMassFractions(initialMassFractions);
    
    // In this example we will not use the policy module (for sake of demonstration of what is happening under the hood)
    //  however, for end users we **strongly** recommend using the policy module to construct engines. It will
    //  ensure that you are not missing important reactions or seed species.

    // 2. Create base network engine (full reaction graph)
    gridfire::GraphEngine baseEngine(comp, NetworkBuildDepth::SecondOrder)

    // 3. Partition network into fast/slow subsets (reduces stiffness)
    gridfire::MultiscalePartitioningEngineView msView(baseEngine);

    // 4. Adaptively cull negligible flux pathways (reduces dimension & stiffness)
    gridfire::AdaptiveEngineView adaptView(msView);

    // 5. Construct implicit solver (handles remaining stiffness)
    gridfire::CVODESolverStrategey solver(adaptView);

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
    std::cout << "Final results are: " << output << std::endl;
}
```

### Callback and Policy Example
Custom callback functions can be registered with any solver. Because it might make sense for each solver to provide
different context to the callback function, you should use the struct `gridfire::solver::<SolverName>::TimestepContext`
as the argument type for the callback function. This struct contains all the information provided by that solver to
the callback function.

```c++
#include "gridfire/gridfire.h" // Unified header for real usage

#include "fourdst/composition/composition.h" // for Composition
#include "fourdst/composition/utils.h" // for buildCompositionFromMassFractions
#include "fourdst/atomic/species.h" // For strongly typed species

#include <iostream>

void callback(const gridfire::solver::CVODESolverStrategy::TimestepContext& context) {
    int H1Index = context.engine.getSpeciesIndex(fourdst::atomic::H_1);
    int He4Index = context.engine.getSpeciesIndex(fourdst::atomic::He_4);

    std::cout << context.t << "," << context.state(H1Index) << "," << context.state(He4Index) << "\n";
}

int main(){
    std::vector<std::string> symbols = {"H-1", "He-4", "C-12"};
    std::vector<double> X = {0.7, 0.29, 0.01};


    const fourdst::composition::Composition composition = fourdst::composition::buildCompositionFromMassFractions(symbols, X);
    gridfire::policy::MainSequencePolicy stellarPolicy(netIn.composition);
    gridfire::engine::DynamicEngine& engine = stellarPolicy.construct();
    
    gridfire::solver::CVODESolverStrategy solver(adaptView);
    solver.set_callback(callback);

    // 6. Prepare input conditions
    gridfire::NetIn input{
        comp,     // composition
        1.5e7,      // temperature [K]
        1.5e2,      // density [g/cm^3]
        1e-12,     // initial timestep [s]
        3e17      // integration end time [s]
    };

    // 7. Execute integration
    gridfire::NetOut output = solver.evaluate(input);
    std::cout << "Final results are: " << output << std::endl;
}
```
>**Note:** If you want to see exactly why each repartitioning stage was triggered in a human readable manner add the flag True to `solver.evaluate` (`solver.evaluate(input, true)`).

>**Note:** A fully detailed list of all available information in the TimestepContext struct is available in the API documentation.

>**Note:** The order of species in the boost state vector (`ctx.state`) is **not guaranteed** to be any particular order run over run. Therefore, in order to reliably extract
> values from it, you **must** use the `getSpeciesIndex` method of the engine to get the index of the species you are interested in (these will always be in the same order).

If you wish to know what is provided by a solver context without investigating the code you can simply do

```c++
void callback(const gridfire::solver::SolverContextBase& context) {
    for (const auto& [parameterName, description] : context.describe()) {
        std::cout << parameterName << ": " << description << "\n";
    }
    std::cout << std::flush();
    exit(0);
}
```

If you set this as the callback (to any solver strategy) it will print out the available parameters and what they
are and then close the code. This is useful when writing new callbacks.


#### Callback Context

Since each solver may provide different context to the callback function, and it may be frustrating to refer to the
documentation every time, we also enforce that all solvers must implement a `descripe_callback_context` method which
returns a vector of tuples<string, string> where the first element is the name of the field and the second is its
datatype. It is on the developer to ensure that this information is accurate.

```c++
...
std::cout << solver.describe_callback_context() << std::endl;
```

## Python
The python bindings intentionally look **very** similar to the C++ code.
Generally all examples can be adapted to python by replacing includes of paths
with imports of modules such that

`#include "gridfire/engine/GraphEngine.h"` becomes `import gridfire.engine.GraphEngine`

All GridFire C++ types have been bound and can be passed around as one would expect.


### Python Example for End Users


The syntax for registration is very similar to C++. There are a few things to note about this more robust example

 1. Note how I use a callback and a log object to store the state of the simulation at each timestep.
 2. If you have tools such as mypy installed you will see that the python bindings are strongly typed. This is
    intentional to help users avoid mistakes when writing code.
```python
from fourdst.composition import Composition
from gridfire.type import NetIn
from gridfire.policy import MainSequencePolicy
from gridfire.solver import CVODESolverStrategy
from enum import Enum
from typing import Dict, Union, SupportsFloat
import json
import dicttoxml

def init_composition() -> Composition:
  Y = [7.0262E-01, 9.7479E-06, 6.8955E-02, 2.5000E-04, 7.8554E-05, 6.0144E-04, 8.1031E-05, 2.1513E-05] # Note these are molar abundances
  S = ["H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"]
  return Composition(S, Y)

def init_netIn(temp: float, rho: float, time: float, comp: Composition) -> NetIn:
  netIn = NetIn()
  netIn.temperature = temp
  netIn.density = rho
  netIn.tMax = time
  netIn.dt0 = 1e-12
  netIn.composition = comp
  return netIn

class StepData(Enum):
  TIME = 0
  DT = 1
  COMP = 2
  CONTRIB = 3


class StepLogger:
  def __init__(self):
    self.num_steps: int = 0
    self.step_data: Dict[int, Dict[StepData, Union[SupportsFloat, Dict[str, SupportsFloat]]]] = {}

  def log_step(self, context):
    engine = context.engine
    self.step_data[self.num_steps] = {}
    self.step_data[self.num_steps][StepData.TIME] = context.t
    self.step_data[self.num_steps][StepData.DT] = context.dt
    comp_data: Dict[str, SupportsFloat] = {}
    for species in engine.getNetworkSpecies():
      sid = engine.getSpeciesIndex(species)
      comp_data[species.name()] = context.state[sid]
    self.step_data[self.num_steps][StepData.COMP] = comp_data
    self.num_steps += 1

  def to_json (self, filename: str):
    serializable_data = {
      stepNum: {
        StepData.TIME.name: step[StepData.TIME],
        StepData.DT.name: step[StepData.DT],
        StepData.COMP.name: step[StepData.COMP],
      }
      for stepNum, step in self.step_data.items()
    }
    with open(filename, 'w') as f:
      json.dump(serializable_data, f, indent=4)

  def to_xml(self, filename: str):
    serializable_data = {
      stepNum: {
        StepData.TIME.name: step[StepData.TIME],
        StepData.DT.name: step[StepData.DT],
        StepData.COMP.name: step[StepData.COMP],
      }
      for stepNum, step in self.step_data.items()
    }
    xml_data = dicttoxml.dicttoxml(serializable_data, custom_root='StepLog', attr_type=False)
    with open(filename, 'wb') as f:
      f.write(xml_data)

def main(temp: float, rho: float, time: float):
  comp = init_composition()
  netIn = init_netIn(temp, rho, time, comp)

  policy = MainSequencePolicy(comp)
  engine = policy.construct()

  solver = CVODESolverStrategy(engine)

  step_logger = StepLogger()
  solver.set_callback(lambda context: step_logger.log_step(context))

  solver.evaluate(netIn, False)
  step_logger.to_xml("log_data.xml")

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description="Simple python example of GridFire usage")
  parser.add_argument("-t", "--temp", type=float, help="Temperature in K", default=1.5e7)
  parser.add_argument("-r", "--rho", type=float, help="Density in g/cm^3", default=1.5e2)
  parser.add_argument("--tMax", type=float, help="Time in s", default=3.1536 * 1e17)
  args = parser.parse_args()
  main(args.temp, args.rho, args.tMax)

```



# Related Projects

GridFire integrates with and builds upon several key 4D-STAR libraries:

- [fourdst](https://github.com/4D-STAR/fourdst): hub module managing versioning
  of `libcomposition`, `libconfig`, `liblogging`, and `libconstants`
- [libcomposition](https://github.com/4D-STAR/libcomposition)
  ([docs](https://4d-star.github.io/libcomposition/)): Composition management
  toolkit.
- [libconfig](https://github.com/4D-STAR/libconfig): Configuration file parsing
  utilities.
- [liblogging](https://github.com/4D-STAR/liblogging): Flexible logging framework.
- [libconstants](https://github.com/4D-STAR/libconstants): Physical constants
