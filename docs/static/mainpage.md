<p align="center">
  <img src="https://github.com/4D-STAR/GridFire/blob/main/assets/logo/GridFire.png?raw=true" width="300" alt="OPAT Core Libraries Logo">
</p>

---
![PyPI - Version](https://img.shields.io/pypi/v/gridfire?style=for-the-badge)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/gridfire?style=for-the-badge)

![GitHub License](https://img.shields.io/github/license/4D-STAR/GridFire?style=for-the-badge)
![ERC](https://img.shields.io/badge/Funded%20by-ERC-blue?style=for-the-badge&logo=europeancommission)

![Dynamic Regex Badge](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fraw.githubusercontent.com%2F4D-STAR%2FGridFire%2Frefs%2Fheads%2Fmain%2Fmeson.build&search=version%3A%20'(.%2B)'%2C&style=for-the-badge&label=Main%20Branch%20Version)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/4D-STAR/GridFire?style=for-the-badge)
![Packaging](https://img.shields.io/github/actions/workflow/status/4D-STAR/GridFire/packaging.yml?branch=main&style=for-the-badge&label=Packaging)


---


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

## Packaging status
Every supported way of installing and linking GridFire is built and tested on each platform wheels are released for,
by the [Packaging workflow](https://github.com/4D-STAR/GridFire/actions/workflows/packaging.yml) (on every push to
`main`). The table shows the latest result on `main`; how each check works, and how to run them locally, is described
in [tests/packaging/ci](tests/packaging/ci/README.md).

| Check | linux-x86_64 | linux-aarch64 | macos-arm64 |
|---|---|---|---|
| wheel builds and repairs as a release does | ![linux-x86_64 wheel](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-x86_64-wheel.json) | ![linux-aarch64 wheel](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-aarch64-wheel.json) | ![macos-arm64 wheel](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/macos-arm64-wheel.json) |
| wheel installs next to the fourdst wheel and imports | ![linux-x86_64 python-install](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-x86_64-python-install.json) | ![linux-aarch64 python-install](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-aarch64-python-install.json) | ![macos-arm64 python-install](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/macos-arm64-python-install.json) |
| C++, C and Fortran programs built with gf-config flags | ![linux-x86_64 gf-config](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-x86_64-gf-config.json) | ![linux-aarch64 gf-config](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-aarch64-gf-config.json) | ![macos-arm64 gf-config](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/macos-arm64-gf-config.json) |
| C++, C and Fortran programs built with the wheel's pkg-config files | ![linux-x86_64 pkg-config](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-x86_64-pkg-config.json) | ![linux-aarch64 pkg-config](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-aarch64-pkg-config.json) | ![macos-arm64 pkg-config](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/macos-arm64-pkg-config.json) |
| concurrent contexts from a GCC/libgomp OpenMP host program | ![linux-x86_64 threads-libgomp](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-x86_64-threads-libgomp.json) | ![linux-aarch64 threads-libgomp](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-aarch64-threads-libgomp.json) | ![macos-arm64 threads-libgomp](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/macos-arm64-threads-libgomp.json) |
| CMake find_package against the wheel (C++, C, Fortran) | ![linux-x86_64 cmake-find](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-x86_64-cmake-find.json) | ![linux-aarch64 cmake-find](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-aarch64-cmake-find.json) | ![macos-arm64 cmake-find](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/macos-arm64-cmake-find.json) |
| CMake FetchContent building GridFire from source (C++, C, Fortran) | ![linux-x86_64 cmake-fetch](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-x86_64-cmake-fetch.json) | ![linux-aarch64 cmake-fetch](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-aarch64-cmake-fetch.json) | ![macos-arm64 cmake-fetch](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/macos-arm64-cmake-fetch.json) |
| Meson subproject from a clean source tree (C++, C, Fortran) | ![linux-x86_64 meson-subproject](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-x86_64-meson-subproject.json) | ![linux-aarch64 meson-subproject](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/linux-aarch64-meson-subproject.json) | ![macos-arm64 meson-subproject](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/4D-STAR/GridFire/ci-status/badges/macos-arm64-meson-subproject.json) |


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
### Currently, known good platforms
The installation script has been tested and found to work on clean
installations of the following platforms:
- MacOS 15.3.2 (Apple Silicon + brew installed)
- Fedora 42.0 (aarch64)
- Ubuntu 24.04 (X86_64)

## Manual Build Instructions

#### Required
- C++ compiler supporting C++23 standard
  - clang >= 17.0.0
  - gcc >= 14.0
- Meson build system (>= 1.5.0)
- Python 3.8 or newer
- CMake 3.20 or newer
- ninja 1.10.0 or newer
- Python packages: `meson-python>=0.15.0`

> **Note:** Windows is not supported at this time and *there are no plans to
> support it in the future*. Windows users are encouraged to use WSL2 or a
> Linux VM.


### Dependency Installation on Common Platforms

- **Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt install pipx pkg-config gcc-14 g++-14
pipx install meson ninja cmake
pipx ensurepath
```
> Note: Ubuntu tends to ship with older versions of compilers, you may need to adjust the specific compiler name
> or add additional repositories if you are running a version of ubuntu older than 24.04 LTS

- **Fedora/CentOS/RHEL:**
```bash
sudo dnf install pkgconf-pkg-config meson ninja-build cmake
```

- **Arch:**
```bash
sudo pacman -Syu pkgconf meson ninja cmake
```

- **macOS (Homebrew):**
```bash
brew update
brew install pkg-config meson ninja cmake gcc
```

> Note: macos >= 15.2.0 ships with a version of clang which works for compilation (and therefore you dont need gcc. I have not tested on versions of macos prior to this).

### Building the C++ Library
```bash
meson setup build
meson compile -C build
```

> Note: On Ubuntu (or any system which does not have a supported version of the compiler as the default) you should use
> ```bash
> CC=gcc-14 CXX=g++-14 meson setup build
> meson compile -C build
> ```

#### Clang vs. GCC
As noted above `clang` tends to compile GridFire much faster than `gcc`. If
your system has both `clang` and `gcc` installed you may force meson to use
clang via environmental variables

> Note: There is an issue I have not resolved on ubuntu 24.04 where clang++-18 fails to find the expected header
> therefore we reccomend the use of gcc-14 and g++-14 on ubuntu 24.04

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
`gcc >= 14.0.0` or `clang >= 17.0.0` work well.


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
- **Solver Module:** Defines the numerical integrators (`PointSolver` for one zone, `GridSolver` for many, both
  CVODE) for solving the stiff ODE systems arising from reaction networks.
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

- **Physics options** live in `gridfire::config::GridFireConfig` (`engine` and `solver` members). Pass one to a
  policy (`MainSequencePolicy(composition, config)`); from C and Fortran use `GFConfig` (see
  [Configuration from C and Fortran](#configuration-from-c-and-fortran)). The engine options are:

  | Option (`config.engine.…`) | Default | Effect |
    |---|---|---|
  | `screening` | `BLENDED` | Plasma screening: `BARE` (none), `WEAK` (Salpeter 1954), `CD09` (Chugunov & DeWitt 2009), or `BLENDED` (Salpeter handing over to CD09 between the `screeningBlend` exponents) |
  | `jacobian` | `ANALYTIC` | How the species Jacobian is computed: `ANALYTIC` or `AUTOMATIC_DIFFERENTIATION` (CppAD) |
  | `useReverseReactions` | `false` | Add reverse (photodisintegration) rates from detailed balance for strong reactions with two reactants and two products. Negligible for hydrogen and helium burning, which is why it is off |
  | `passiveSpecies` | none | Species carried in the composition without reactions: they count towards Y_e and screening but never change |
  | `energyChannels.splitHydrogenPPCNO` | `true` | Report hydrogen burning as separate pp and CNO energy channels |

### Rate multipliers
Any reaction's rate can be scaled by a constant factor, for example to sample rate uncertainties in Monte Carlo
runs. A factor applies to everything computed afterwards: abundances, energy generation, neutrino losses, both
Jacobian methods, the eps derivatives and the per-reaction contributions. When reverse reactions are on, the reverse
rate is scaled with its forward rate, so detailed balance holds. The network is not rebuilt, so changing factors
between runs is cheap.

Reactions are named by the ids GridFire reports (`getNetworkReactions()` in C++ and Python, `gf_get_reaction_id` in
C, `gff_get_reaction_ids` in Fortran), for example `p(p,e+)d` or `n14(p,g)o15`. A factor of 1 removes a multiplier
and 0 switches the reaction off. An unknown reaction id, or a negative or non-finite factor, is an error, so a host
code can check its list of perturbed rates at start-up.

| | Set | Read | Clear |
|---|---|---|---|
| C++ / Python (`GraphEngine`) | `setRateMultiplier(id, factor)` | `getRateMultiplier(id)`, `getRateMultipliers()` | `clearRateMultipliers()` |
| C | `gf_set_rate_multiplier(ctx, id, factor)` | `gf_get_rate_multiplier(ctx, id, &factor)` | `gf_clear_rate_multipliers(ctx)` |
| Fortran | `net%gff_set_rate_multiplier(id, factor[, ierr])` | `net%gff_get_rate_multiplier(id, factor[, ierr])` | `net%gff_clear_rate_multipliers([ierr])` |

In C++ and Python the `GraphEngine` is the bottom of a policy's engine stack:

```python
policy = MainSequencePolicy(comp)
net = policy.construct()
graph = policy.get_engine_stack()[0]
graph.setRateMultiplier("p(p,e+)d", 1.1)      # +10% on p + p
graph.setRateMultiplier("n14(p,g)o15", 0.9)   # -10% on N-14(p,g)
print(graph.getRateMultipliers())             # [('n14(p,g)o15', 0.9), ('p(p,e+)d', 1.1)]
```

Do not change multipliers while another thread is evaluating with the same engine.

### Provenance
`gf_get_provenance` (C), `gff_get_provenance` (Fortran) and `gridfire::describe_provenance` (C++) describe what a
calculation uses, one `key: value` line per item, for a host code to write next to its results:

```
gridfire_version: 1.2.0
policy: MAIN_SEQUENCE_POLICY
network: 98 species, 350 reactions, hash 1a113533eac1d97d
strong_rates: REACLIB (pynucastro snapshot reaclib_pynucastro_latest_Z26), data hash 5ab8a8f6f54ff09c
weak_rates: Ravlic, Giraud, Paar & Zegers (2024), arXiv:2412.00650, data hash 60af6d0fc70e66b2
rate_multipliers: n14(p,g)o15=0.9, p(p,e+)d=1.1
screening: BLENDED (Salpeter 1954 -> Chugunov & DeWitt 2009), handover for weak-screening exponent 0.5 to 1.5
jacobian: ANALYTIC
reverse_reactions: off
passive_species: none
energy_channels: hydrogen split into pp and CNO
solver: CVODE, abs_tol 1e-08, rel_tol 1e-05, clamp_negative_iterates false, max_step_failure_restarts 3
```

Two runs that print the same description used the same code version, rate data, network and physics options.
`gf_get_version()` / `gff_get_version()` give the version alone.

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

GridFire integrates a network with a solver that works on any `DynamicEngine`. The solver itself is stateless and
can be shared; the state of one calculation (CVODE memory, step history, the engine's scratch space) lives in a
separate *context*, so one engine and solver can serve many zones or threads, each with its own context.

### Solvers
- `gridfire::solver::PointSolver`: one zone. `evaluate(ctx, netIn)` integrates to `netIn.tMax` and returns a `NetOut`;
  `evaluateRates(ctx, netIn)` returns the instantaneous rates (`InstantaneousRates`) without integrating.
- `gridfire::solver::GridSolver`: many zones at once (spread over OpenMP threads when available), built on a
  `PointSolver`.

The context for a `PointSolver` is a `PointSolverContext`, made from the engine's state blob (the `scratch_blob` a
policy's `construct()` returns, or `engine.constructStateBlob(nullptr)`).

### NetIn and NetOut
GridFire solvers use a unified input and output type for their public interface.

#### NetIn
A `NetIn` struct contains
- The composition to start the timestep at, as molar abundances (`NetIn::composition`)
- The temperature in Kelvin (`NetIn::temperature`)
- The density in g/cm^3 (`NetIn::density`)
- The time to evolve the network for, in seconds (`NetIn::tMax`)
- The initial timestep in seconds (`NetIn::dt0`)
- The initial energy in the system in erg/g (`NetIn::energy`)

>**Note:** The composition must be a `fourdst::composition::Composition`
>object, from the `fourdst/composition/composition.h` header. `fourdst` is installed
>automatically with GridFire. In Python it is `fourdst.composition.Composition`.

#### NetOut
A `NetOut` struct contains
- The final composition (`NetOut::composition`) and the number of steps taken (`NetOut::num_steps`)
- The specific energy released over the step (`NetOut::energy`) [erg/g]
- The neutrino energy lost and neutrinos emitted over the step (`NetOut::specific_neutrino_energy_loss`,
  `NetOut::specific_neutrino_flux`)
- At the final state: the energy generation rate net of neutrinos (`NetOut::eps_nuc`) [erg/g/s], its derivatives
  (`dEps_dT`, `dEps_dRho`, `dlnEps_dlnT`, `dlnEps_dlnRho`), the neutrino loss rate
  (`neutrino_energy_loss_rate`) and `eps_nuc` split by energy channel (`energy_channel_names`, `eps_by_channel`)
- For multi-zone solvers, whether the zone converged (`converged`, `time_reached`, `failure_message`)

### PointSolver

We use the CVODE module from [SUNDIALS](https://computing.llnl.gov/projects/sundials/cvode) as our numerical solver,
specifically its BDF linear multistep method with adaptive timestepping and a dense linear solver. The tolerances
(`solver.pointSolver.absTol`, `relTol`) and the step-failure handling are set in the configuration.

A callback can be registered, in Python or C++, on the solver context (`PointSolverContext::callback`). It is called
after every accepted step with a `const PointSolverTimestepContext&`, which holds the time, step size, state vector,
engine and more. Callbacks are for inspection: they cannot change the state of the network. Physics changes belong
in an engine or engine view.

## Python Extensibility
Through the Python bindings, users can subclass engine view classes directly in
Python, override methods like `evaluate` or `generateStoichiometryMatrix`, and
pass instances back into C++ solvers. This enables rapid prototyping of custom
view strategies without touching C++ sources.


# Usage Examples

## C++
All of GridFire is available through the unified header `gridfire/gridfire.h`. See
[Building a program against GridFire](#building-a-program-against-gridfire) for compiler flags. For example, with the
pip wheel installed:

```bash
c++ -std=c++23 example.cpp $(gf-config --cflags --libs) -o example
```

### Composition Initialization
```c++
#include "fourdst/composition/composition.h"
#include "fourdst/composition/utils/utils.h" // buildCompositionFromMassFractions

#include <iostream>
#include <string>
#include <vector>

int main() {
    // From molar abundances (what GridFire works in) ...
    const fourdst::composition::Composition fromY(
        std::vector<std::string>{"H-1", "He-4", "C-12"}, {0.7, 0.0725, 8.3e-4});

    // ... or from mass fractions, which must sum to 1
    const fourdst::composition::Composition fromX = fourdst::composition::buildCompositionFromMassFractions(
        std::vector<std::string>{"H-1", "He-4", "C-12"}, {0.7, 0.29, 0.01});

    std::cout << fromY << "\n" << fromX << "\n";
    std::cout << "X(He-4) = " << fromX.getMassFraction("He-4") << ", Y(He-4) = " << fromX.getMolarAbundance("He-4") << "\n";
}
```

### Common Workflow Example
Most programs let a policy build the network, then burn with a `PointSolver`:

```c++
#include "gridfire/gridfire.h"
#include "fourdst/composition/composition.h"

#include <iostream>
#include <string>
#include <vector>

int main() {
    using namespace gridfire;

    // 1. Initial composition (molar abundances)
    const fourdst::composition::Composition comp(
        std::vector<std::string>{"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"},
        {0.7026166, 9.7479e-06, 0.0689551, 0.00025, 7.8554e-05, 0.00060144, 8.103e-05, 2.1513e-05});

    // 2. Let a policy build the network: it guarantees the seed species and key reactions are present
    config::GridFireConfig config;                  // defaults; change physics options here
    policy::MainSequencePolicy policy(comp, config);
    const policy::ConstructionResults net = policy.construct();

    // 3. A solver for the engine, and a context holding one calculation's state (CVODE memory, scratch space)
    const solver::PointSolver solver(net.engine);
    solver::PointSolverContext ctx(*net.scratch_blob);
    ctx.stdout_logging = false;                     // no per-step log

    // 4. Burn for 1 Myr at solar-core conditions
    NetIn input;
    input.composition = comp;
    input.temperature = 1.5e7;  // K
    input.density = 150.0;      // g/cm^3
    input.tMax = 3.15e13;       // s
    input.dt0 = 1e-12;          // s, first step (CVODE grows it quickly)

    const NetOut output = solver.evaluate(ctx, input);
    std::cout << output.num_steps << " steps; X(H-1): " << comp.getMassFraction("H-1") << " -> "
              << output.composition.getMassFraction("H-1") << "\n"
              << "eps_nuc = " << output.eps_nuc << " erg/g/s, energy released = " << output.energy << " erg/g\n";

    // 5. Instantaneous rates at a state, without integrating
    const InstantaneousRates rates = solver.evaluateRates(ctx, input);
    std::cout << "eps_nuc now = " << rates.eps_nuc << " erg/g/s, d eps/dT = " << rates.dEps_dT << "\n";
}
```

### GraphEngine Initialization
A `GraphEngine` can also be built directly, from a seed composition or an explicit `reaction::ReactionSet`. Unlike a
policy, this does not check that the network contains the species and reactions a burning stage needs.

```c++
#include "gridfire/gridfire.h"
#include "fourdst/composition/composition.h"

#include <iostream>
#include <string>
#include <vector>

int main() {
    using namespace gridfire;
    const fourdst::composition::Composition comp(
        std::vector<std::string>{"H-1", "He-4", "C-12", "N-14", "O-16"}, {0.7, 0.07, 2.5e-4, 7.9e-5, 6.0e-4});

    // Physics options for the engine (see "GraphEngine Configuration Options")
    config::GridFireConfig config;
    config.engine.screening = screening::ScreeningType::BLENDED;
    config.engine.useReverseReactions = false;

    // Grow the network from the seed composition: SecondOrder follows reactions two steps out
    const partition::RauscherThielemannPartitionFunction partitionFunction; // used for reverse rates
    const engine::GraphEngine engine(comp, partitionFunction,
                                     engine::NetworkBuildDepth::SecondOrder,
                                     engine::NetworkConstructionFlags::DEFAULT, config);

    // Engines keep per-calculation state in a StateBlob, so one engine can serve many threads or zones
    const auto ctx = engine.constructStateBlob(nullptr);
    std::cout << engine.getNetworkSpecies(*ctx).size() << " species, "
              << engine.getNetworkReactions(*ctx).size() << " reactions\n";
}
```

### Callback Example
```c++
#include "gridfire/gridfire.h"
#include "fourdst/composition/composition.h"
#include "fourdst/atomic/species.h" // strongly typed species, e.g. fourdst::atomic::H_1

#include <nvector/nvector_serial.h> // NV_Ith_S, to read the CVODE state vector

#include <iostream>
#include <string>
#include <vector>

void callback(const gridfire::solver::PointSolverTimestepContext& context) {
    // The state vector holds the network species in the engine's order (then the accumulated energy), so look
    // species up by index rather than assuming an order
    const size_t iH1 = context.engine.getSpeciesIndex(context.state_ctx, fourdst::atomic::H_1);
    const size_t iHe4 = context.engine.getSpeciesIndex(context.state_ctx, fourdst::atomic::He_4);
    if (context.num_steps % 50 == 0) {
        std::cout << context.t << "," << NV_Ith_S(context.state, iH1) << "," << NV_Ith_S(context.state, iHe4) << "\n";
    }
}

int main() {
    using namespace gridfire;
    const fourdst::composition::Composition comp(
        std::vector<std::string>{"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"},
        {0.7026166, 9.7479e-06, 0.0689551, 0.00025, 7.8554e-05, 0.00060144, 8.103e-05, 2.1513e-05});

    policy::MainSequencePolicy policy(comp);
    const policy::ConstructionResults net = policy.construct();
    const solver::PointSolver solver(net.engine);

    solver::PointSolverContext ctx(*net.scratch_blob);
    ctx.stdout_logging = false;
    ctx.callback = callback;      // called after every accepted step

    NetIn input;
    input.composition = comp;
    input.temperature = 1.5e7;
    input.density = 150.0;
    input.tMax = 3.15e13;
    input.dt0 = 1e-12;

    std::cout << "t,Y(H-1),Y(He-4)\n";
    const NetOut output = solver.evaluate(ctx, input);
    std::cout << output << "\n";
}
```

>**Note:** The order of species in the state vector (`context.state`) is the engine's, not the order of the input
>composition. Use the engine's `getSpeciesIndex` (or `context.networkSpecies`) to find a species.

To see everything a timestep context provides, register a callback that prints its `describe()` and exits:

```c++
void describe_callback(const gridfire::solver::PointSolverTimestepContext& context) {
    for (const auto& [name, description] : context.describe()) {
        std::cout << name << ": " << description << "\n";
    }
    std::exit(0);
}
```

## Python
The python bindings intentionally look **very** similar to the C++ code: namespaces become modules, so
`gridfire::policy::MainSequencePolicy` is `gridfire.policy.MainSequencePolicy`, `gridfire::solver::PointSolver` is
`gridfire.solver.PointSolver`, and `gridfire::NetIn` is `gridfire.type.NetIn`. Compositions come from
`fourdst.composition`.

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
from gridfire.solver import PointSolver, PointSolverContext
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
    comp = context.composition
    for species in engine.getNetworkSpecies(context.state_ctx):
      comp_data[species.name()] = comp.getMolarAbundance(species)
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
  construct = policy.construct()

  engine = construct.engine
  context = PointSolverContext(construct.scratch_blob)

  solver = PointSolver(engine)

  logger = StepLogger()

  context.callback = lambda ctx: logger.log_step(ctx)
  solver.evaluate(context, netIn)
  logger.to_xml("log_data.xml")

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description="Simple python example of GridFire usage")
  parser.add_argument("-t", "--temp", type=float, help="Temperature in K", default=1.5e7)
  parser.add_argument("-r", "--rho", type=float, help="Density in g/cm^3", default=1.5e2)
  parser.add_argument("--tMax", type=float, help="Time in s", default=3.1536 * 1e17)
  args = parser.parse_args()
  main(args.temp, args.rho, args.tMax)
```


# External Usage
C++ does not have a stable ABI nor does it make any strong guarantees about stl container layouts between compiler versions.
Therefore, GridFire includes a set of stable C bindings which can be used to interface with a limited subset of GridFire functionality
from other languages.

> **Note:** These bindings are not intended to allow GridFire to be extended from other languages; rather, they are intended to allow GridFire to be used as a
> black-box library from other languages.

> **Note:** One assumption for external usage is that the ordering of the species list will not change. That is to say that whatever order the array
> used to register the species is will be assumed to always be the order used when passing abundance arrays to and from GridFire.

> **Note:** Because the C API does not pass the general Composition object a `mass_lost`
> output parameter has been added to the evolve calls, this tracks the total mass in species which have not been registered with the C API GridFire by the caller

> **Note:** The C and Fortran APIs take and return **molar abundances** (Y). `gf_mass_fractions_to_molar` and
> `gf_molar_to_mass_fractions` (`gff_mass_fractions_to_molar` / `gff_molar_to_mass_fractions`) convert, in the
> registered species order.

## Building a program against GridFire
The pip wheel ships everything a C, C++ or Fortran program needs: headers, `libgridfire`, `libgridfire_extern` (the
C API), `gridfire_mod.f90` (the Fortran module source, compiled with your own Fortran compiler so its `.mod` file
matches), pkg-config files and a CMake package.

**`gf-config`**, installed with the wheel, prints compiler and linker flags:

```bash
c++ -std=c++23 prog.cpp $(gf-config --cflags --libs) -o prog                               # C++
cc prog.c $(gf-config --component c --cflags --libs) -o prog                               # C API
gfortran $(gf-config --fortran-source) prog.f90 $(gf-config --component fortran --cflags --libs) -o prog  # Fortran
```

**pkg-config:** the packages are `gridfire`, `gridfire_extern` and `gridfire_fortran`. The Fortran package's
`fortran_source` variable is the path of `gridfire_mod.f90`.

```bash
export PKG_CONFIG_PATH="$(gf-config --pkg-config-path)"
cc prog.c $(pkg-config --cflags --libs gridfire_extern) -o prog
```

**CMake, with an installed GridFire** (the wheel, or `meson install`): the imported targets are
`GridFire::gridfire` (C++), `GridFire::gridfire_extern` (C) and `GridFire::gridfire_fortran` (Fortran; it compiles
`gridfire_mod.f90` with your project's compiler). Configure with
`-DCMAKE_PREFIX_PATH="$(gf-config --cmake-prefix)"` for the wheel.

```cmake
find_package(GridFire 1.2 REQUIRED COMPONENTS CXX C Fortran)
target_link_libraries(my_code PRIVATE GridFire::gridfire_fortran)
```

**CMake, building GridFire from source:** FetchContent uses an installed GridFire when it finds one and otherwise
fetches it and builds it with Meson (which, with ninja, must be installed) while CMake configures:

```cmake
include(FetchContent)
FetchContent_Declare(GridFire
    GIT_REPOSITORY https://github.com/4D-STAR/GridFire.git
    GIT_TAG        v1.2.0
    GIT_SHALLOW    TRUE
    FIND_PACKAGE_ARGS 1.2 COMPONENTS CXX C Fortran)
FetchContent_MakeAvailable(GridFire)
```

The cache variables `GRIDFIRE_BUILDTYPE` and `GRIDFIRE_MESON_OPTIONS` pass options through to that Meson build.

**Meson subproject:** put [`tests/packaging/meson_consumer/gridfire.wrap`](tests/packaging/meson_consumer/gridfire.wrap)
in your `subprojects/` directory, then

```meson
gridfire_dep = dependency('gridfire', version: '>=1.2.0', default_options: ['build_fortran=true', 'build_tests=false'])
gridfire_extern_dep = dependency('gridfire_extern', version: '>=1.2.0')
gridfire_fortran_dep = dependency('gridfire_fortran', version: '>=1.2.0')
```

Each route is checked on every supported platform; see [Packaging status](#packaging-status).

## Thread safety
In OpenMP builds, which include the published wheels:

- Build and free contexts from serial code, while no other GridFire call is running: `gf_init` through
  `gf_construct_solver_from_engine`, and `gf_free`. Make the first `gf_init` call outside any parallel region.
- Once built, **separate** contexts may evolve and evaluate concurrently from any threads, including OpenMP threads of
  a GCC/libgomp host program such as a `gfortran -fopenmp` code.
- A single context must never be used by two threads at once.
- `MULTI_ZONE` contexts spread their zones over GridFire's own OpenMP threads; call them from serial code.

## C API Overview
In general when using the C API the workflow is to

1. Create a context with `gf_init(SINGLE_ZONE)` (or `MULTI_ZONE`, then `gf_set_num_zones`).
2. Optionally set the configuration (`gf_config_defaults`, then `gf_set_config`, or `gf_load_config` from TOML).
3. Register your species, build the engine with `gf_construct_engine_from_policy`, then the solver with
   `gf_construct_solver_from_engine`.
4. Call `gf_evolve` to burn over a time step, `gf_evaluate` for instantaneous rates without integrating, or
   `gf_evaluate_reactions` for per-reaction contributions.
5. Check every return code: `GF_SUCCESS` (0) is success (`gf_register_species` returns `FDSSE_SUCCESS`, 1), and
   `gf_get_last_error_message(ctx)` describes a failure.
6. Call `gf_free` to release the context.

### Configuration from C and Fortran
`struct GFConfig` (`type(GFConfig)` in Fortran) holds the solver tolerances and engine options. Fill it with the
defaults first, change what you need, and set it **before** `gf_construct_engine_from_policy`:

```C
struct GFConfig config;
gf_config_defaults(&config);
config.engine_screening = GF_SCREENING_CD09;
config.engine_jacobian = GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION;
config.engine_use_reverse_reactions = 1;
config.solver_rel_tol = 1e-8;
gf_set_config(ctx, &config);
```

`gf_save_config` and `gf_load_config` write and read the same settings as TOML.

### C Example

```C
#include "gridfire/extern/gridfire_extern.h"
#include <stdio.h>

#define NUM_SPECIES 8

/* Report a failed call and stop (gf_register_species reports success as FDSSE_SUCCESS) */
#define CHECK(call, ctx, what)                                                                          \
    do {                                                                                                \
        const int ret_ = (call);                                                                        \
        if (ret_ != GF_SUCCESS && ret_ != FDSSE_SUCCESS) {                                              \
            printf("%s failed: %s [%s]\n", what, gf_get_last_error_message(ctx), gf_error_code_to_string(ret_)); \
            gf_free(SINGLE_ZONE, ctx);                                                                  \
            return 1;                                                                                   \
        }                                                                                               \
    } while (0)

int main(void) {
    const char* species_names[NUM_SPECIES] = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
    /* Molar abundances */
    const double Y[NUM_SPECIES] = {0.7026166, 9.7479e-06, 0.0689551, 0.00025, 7.8554e-05, 0.00060144, 8.103e-05, 2.1513e-05};
    const double T = 1.5e7;    /* K */
    const double rho = 150.0;  /* g/cm^3 */

    void* ctx = gf_init(SINGLE_ZONE);
    CHECK(gf_register_species(ctx, NUM_SPECIES, species_names), ctx, "gf_register_species");
    CHECK(gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", Y, NUM_SPECIES), ctx, "engine");
    CHECK(gf_construct_solver_from_engine(ctx), ctx, "solver");

    /* Perturb two rates, e.g. for one Monte Carlo sample */
    CHECK(gf_set_rate_multiplier(ctx, "p(p,e+)d", 1.1), ctx, "gf_set_rate_multiplier");
    CHECK(gf_set_rate_multiplier(ctx, "n14(p,g)o15", 0.9), ctx, "gf_set_rate_multiplier");

    /* Outputs: the caller provides arrays of the right size */
    size_t num_channels = 0;
    CHECK(gf_get_num_energy_channels(ctx, &num_channels), ctx, "gf_get_num_energy_channels");
    double Y_out[NUM_SPECIES], eps_by_channel[16];
    double energy, dEps_dT, dEps_dRho, nu_loss, nu_flux, mass_lost, eps_nuc, dlnEps_dlnT, dlnEps_dlnRho, nu_rate;
    if (num_channels > 16) return 1;

    CHECK(gf_evolve(SINGLE_ZONE, ctx, Y, NUM_SPECIES, &T, &rho,
                    3.15e13,  /* burn for 1 Myr (s) */
                    1e-12,    /* initial step (s) */
                    Y_out, &energy, &dEps_dT, &dEps_dRho, &nu_loss, &nu_flux, &mass_lost,
                    &eps_nuc, &dlnEps_dlnT, &dlnEps_dlnRho, &nu_rate, eps_by_channel),
          ctx, "gf_evolve");

    for (size_t i = 0; i < NUM_SPECIES; i++) {
        printf("%-6s %.6e -> %.6e\n", species_names[i], Y[i], Y_out[i]);
    }
    printf("eps_nuc = %.6e erg/g/s, energy released = %.6e erg/g\n", eps_nuc, energy);

    const char* provenance = NULL;
    CHECK(gf_get_provenance(ctx, &provenance), ctx, "gf_get_provenance");
    printf("%s", provenance);

    gf_free(SINGLE_ZONE, ctx);
    return 0;
}
```

## Fortran API Overview

GridFire makes use of the stable C API and Fortran 2003's `iso_c_binding` to provide a Fortran interface for legacy
code. The Fortran interface mirrors the C API: every C function has a type-bound procedure with the `gff_` prefix
(GridFire Fortran) on a `type(GridFire)` object.

1. `net%gff_init(SINGLE_ZONE)` (or `gff_init(MULTI_ZONE, zones)`): create the context.
2. `net%gff_set_config(config)`: optional; `config = gff_default_config()`, then change fields of the `GFConfig`.
3. `net%gff_register_species(names)`, `net%gff_setup_policy("MAIN_SEQUENCE_POLICY", Y)`, `net%gff_setup_solver()`.
4. `net%gff_evolve(...)`, `net%gff_evaluate(...)`, `net%gff_evaluate_reactions(...)`: scalar arguments for one zone,
   arrays for `MULTI_ZONE`.
5. `net%gff_set_rate_multiplier`, `net%gff_get_provenance()`, `gff_get_version()`: see
   [Rate multipliers](#rate-multipliers) and [Provenance](#provenance).
6. `net%gff_get_last_error()`: the message for the last failure.
7. `net%gff_free()`.

Calls that can fail take an `ierr` argument. Where it is optional (the configuration and rate-multiplier calls), a
failure stops the program when `ierr` is absent and is returned when it is present.

The module is shipped as source, `gridfire_mod.f90`, because `.mod` files are specific to a compiler and its version:
compile it with the compiler that builds your program and link against `libgridfire_extern`. `gf-config`,
pkg-config, CMake and Meson all do this for you; see [Building a program against GridFire](#building-a-program-against-gridfire).

### Fortran Example

```fortran
program main
  use iso_c_binding
  use gridfire_mod
  implicit none

  ! Names are blank-padded; the module trims them
  character(len=5), dimension(8) :: species_names = [ &
          "H-1  ", "He-3 ", "He-4 ", "C-12 ", "N-14 ", "O-16 ", "Ne-20", "Mg-24" ]
  ! Molar abundances (Y = X / A), solar-ish
  real(c_double), dimension(8) :: Y_in = [ &
          0.7026166d0, 9.7479d-06, 0.0689551d0, 0.00025d0, 7.8554d-05, 0.00060144d0, 8.103d-05, 2.1513d-05 ]
  real(c_double), dimension(8) :: Y_out
  real(c_double), allocatable :: eps_by_channel(:)
  real(c_double) :: energy, dedt, dedrho, nu_loss, nu_flux, mass_lost, eps_nuc, dlneps_dlnt, dlneps_dlnrho, nu_rate
  real(c_double) :: T = 1.5d7       ! K
  real(c_double) :: rho = 150.0d0   ! g/cm^3
  real(c_double) :: tMax = 3.15d13  ! 1 Myr, in s
  type(GridFire) :: net
  type(GFConfig) :: config
  integer :: i, ierr

  call net%gff_init(SINGLE_ZONE)

  config = gff_default_config()
  config%engine_screening = GF_SCREENING_BLENDED
  call net%gff_set_config(config)

  call net%gff_register_species(species_names)
  call net%gff_setup_policy("MAIN_SEQUENCE_POLICY", Y_in)
  call net%gff_setup_solver()

  ! One Monte Carlo sample: perturb two rates. With ierr present an unknown id is reported, not fatal
  call net%gff_set_rate_multiplier("p(p,e+)d", 1.1d0, ierr)
  if (ierr /= GF_SUCCESS) print *, "p(p,e+)d is not in the network"
  call net%gff_set_rate_multiplier("n14(p,g)o15", 0.9d0)

  allocate(eps_by_channel(net%gff_get_num_energy_channels()))
  call net%gff_evolve(Y_in, T, rho, tMax, 1d-12, Y_out, energy, dedt, dedrho, nu_loss, nu_flux, mass_lost, &
          eps_nuc, dlneps_dlnt, dlneps_dlnrho, nu_rate, eps_by_channel, ierr)
  if (ierr /= GF_SUCCESS) then
    print *, "Evolution failed: ", net%gff_get_last_error()
    call net%gff_free()
    error stop 1
  end if

  print '(a, es13.6, a)', "eps_nuc: ", eps_nuc, " erg/g/s"
  do i = 1, size(species_names)
    print '(a6, es13.6, " -> ", es13.6)', species_names(i), Y_in(i), Y_out(i)
  end do
  print '(a)', net%gff_get_provenance()

  call net%gff_free()
end program main
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
- [libplugin](https://github.com/4D-STAR/libplugin): Dynamically loadable plugin
  framework.
