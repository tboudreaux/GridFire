
![GridFire Logo](../../assets/logo/GridFire.png)

# GridFire: A Nuclear Reaction Network Library

## Overview

GridFire is a C++ library for simulating nuclear reaction networks. It is designed to be flexible, performant, and easy to integrate into larger astrophysical simulations. GridFire provides a modular framework for defining, manipulating, and solving complex reaction networks, with a focus on providing different "engines" to suit various computational needs.

This documentation provides a comprehensive overview of the GridFire library, its architecture, and its components.

## Features

* **Modular Design:** The library is organized into distinct components for reactions, engines, and solvers, allowing for easy extension and customization.
* **Multiple Engines:** GridFire offers several reaction network engines, including:
    * A `GraphEngine` for solving the full reaction network.
    * An `AdaptiveEngine` that dynamically adjusts the network size at runtime for improved performance.
    * A simplified `Approx8Engine` for testing and for scenarios where a full network is not required.
* **REACLIB Integration:** GridFire can parse and utilize data from the REACLIB database, allowing for the use of standard, up-to-date reaction rates.
* **High Performance:** The library is designed for performance, using modern C++ features and high-performance libraries like Eigen, Boost, and CppAD.

## Directory Structure

The GridFire project is organized into the following main directories:

* `src/`: Contains the source code for the GridFire library.
    * `include/gridfire/`: Public header files for the library.
        * `reaction/`: Header files related to defining and managing nuclear reactions.
        * `engine/`: Header files for the various reaction network engines.
        * `solver/`: Header files for the numerical solvers.
    * `lib/`: Source code for the implementation of the library components.
* `subprojects/`: Contains external libraries and dependencies.
* `tests/`: Unit tests for the GridFire library.
* `docs/`: Directory for documentation, including this Doxygen homepage.

## Core Components

The GridFire library is built around three main components: **Reactions**, **Engines**, and **Solvers**.

### Reactions

The `reaction` component is responsible for defining and managing nuclear reactions. The key classes are:

* `gridfire::reaction::Reaction`: Represents a single nuclear reaction, including its reactants, products, and rate coefficients.
* `gridfire::reaction::ReactionSet`: A collection of `Reaction` objects, representing a full reaction network.
* `gridfire::reaction::LogicalReaction`: An abstraction that can represent a single reaction or a combination of multiple reactions (e.g., a forward and reverse reaction).
* `gridfire::reaction::LogicalReactionSet`: A collection of `LogicalReaction` objects.

GridFire can load reaction data from the REACLIB database via the `build_reaclib_nuclear_network` function, which is implemented in `src/network/lib/reaction/reaclib.cpp`.

### Engines

Engines are responsible for the low-level details of solving the reaction network equations. GridFire provides several engines, each with its own trade-offs in terms of accuracy and performance:

* **`gridfire::GraphEngine`**: A full-network solver that represents the reaction network as a graph. It uses the CppAD library for automatic differentiation to compute the Jacobian matrix, which is essential for implicit solvers.
* **`gridfire::AdaptiveEngine`**: A dynamic engine that adapts the size of the reaction network at runtime. It can add or remove reactions and species based on their importance, which can significantly improve performance in simulations with changing conditions.
* **`gridfire::Approx8Engine`**: A simplified, 8-isotope network for hydrogen and helium burning. This engine is useful for testing and for simulations where a full network is not necessary. It includes the pp-chain and the CNO cycle.

### Solvers

The `solver` component provides the numerical algorithms for solving the system of ordinary differential equations (ODEs) that describe the evolution of the species abundances. GridFire uses the `boost::numeric::odeint` library for its ODE solvers.

The main solver class is `gridfire::solver::QSENetworkSolver`, which implements a Quasi-Steady-State (QSE) solver.

## Dependencies

GridFire relies on the following external libraries:

* **[Eigen](https://eigen.tuxfamily.org/)**: A C++ template library for linear algebra.
* **[Boost](https://www.boost.org/)**: A collection of peer-reviewed, high-quality C++ libraries. GridFire uses Boost for sparse matrices and numerical integration.
* **[CppAD](https://cppad.readthedocs.io/en/latest/)**: A C++ template library for automatic differentiation.
* **[Quill](https://github.com/odygrd/quill)**: A high-performance, asynchronous logging library.
* **[yaml-cpp](https://github.com/jbeder/yaml-cpp)**: A YAML parser and emitter for C++.
* **[GoogleTest](https://github.com/google/googletest)**: A unit testing framework for C++.
* **[xxHash](https://github.com/Cyan4973/xxHash)**: An extremely fast non-cryptographic hash algorithm.

Aside from boost, all dependencies are automatically downloaded and built by the Meson build system. Boost must be installed on your system beforehand.

## Building the Project

GridFire uses the [Meson](https://mesonbuild.com/) build system. To build the project, you will need to have Meson and a C++ compiler (like g++ or clang++) installed.

You can configure the build using the options listed in `meson_options.txt`.

```bash
meson setup builddir --prefix=/path/to/install --buildtype=release
meson compile -C builddir
meson install -C builddir
```

## How to Use

To use the GridFire library in your own project, you will need to link against the compiled library and include the necessary header files.

Here is a simple example of how to use the library:

```cpp
#include "gridfire/network.h"
#include "gridfire/solver/solver.h"

int main() {
    // 1. Create a Composition object with the initial abundances
    fourdst::composition::Composition comp;
    // ... set initial abundances ...

    // 2. Create a reaction network engine
    gridfire::GraphEngine engine(comp);

    // 3. Create a solver
    gridfire::solver::QSENetworkSolver solver(engine);

    // 4. Set the thermodynamic conditions
    gridfire::NetIn netIn;
    netIn.temperature = 1.0e8; // K
    netIn.density = 1.0e4; // g/cm^3

    // 5. Solve the network
    gridfire::NetOut netOut = solver.evaluate(netIn);

    return 0;
}
```
Linking can be done using the pkg-config tool:

```bash
pkg-config --cflags --libs gridfire
```

Note that you will also need to tell the compiler where to find boost headers.

A more detailed example of how to use an AdaptiveEngine can be found below
```cpp
#include "gridfire/network.h"
#include "gridfire/solver/solver.h"

int main() {
    // 1. Create a Composition object with the initial abundances
    fourdst::composition::Composition comp;
    // ... set initial abundances ...

    gridfire::GraphEngine baseEngine(comp);
    gridfire::AdaptiveEngineView adaptiveEngine(comp);

    // 3. Create a solver
    gridfire::solver::QSENetworkSolver solver(adaptiveEngine);

    // 4. Set the thermodynamic conditions
    gridfire::NetIn netIn;
    netIn.temperature = 1.0e8; // K
    netIn.density = 1.0e4; // g/cm^3

    // 5. Solve the network
    gridfire::NetOut netOut = solver.evaluate(netIn);

    return 0;
}
```

Note how the adaptive engine is a view of the base engine. This allows the adaptive engine to dynamically adjust the
network size at runtime based on the reactions that are active.