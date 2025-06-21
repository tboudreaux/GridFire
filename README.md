<p align="center">
  <img src="assets/logo/GridFire.png" width="300" alt="OPAT Core Libraries Logo">
</p>


---

GridFire is a C++ library designed to preform general nuclear network evolution using the reaclib library. It is a
part of the larger SERiF project which is itself a part of the 4D-STAR collaboration. 

> Note that GridFire is still very early in development and is not ready for scientific use. The API is not stable and
> the reults are not yet validated. We are actively working on improving the library and once it is tested and validated
> we will release a stable version.

## Design

GridFire uses a "graph-first" design to represent the nuclear network. Specifically, internally the network is
represented as a directed hypergraph where nodes are nuclides and edges are reactions. This design allows for very straightforward
and efficient evolution of the network topology. 

## Current Features
- **Reaclib Support**: GridFire uses reaclib coefficients statically compiled into the binary to evaluate reation rates. We bundle a script which can be used to generate the header files where these are stored. 
- **Auto Differentiation**: GridFire uses [CppAD](https://github.com/coin-or/CppAD) to generate analytic Jacobians for stiff network evaluation.
- **Dynamic Network Topology**: GridFire supports dynamic network topology, allowing for the addition and removal of nuclides and reactions during runtime.
- **Dynamic Stiff Detection**: GridFire uses a heuristic stiff detection algorithm to determine the stiffness of the network at runtime.

## Planned Features
### High Priority
- **Reverse Reactions**: Currently, GridFire only supports forward reactions. Very soon we will add support for reverse reactions.
- **Weak Reactions**: GridFire will eventually support weak reactions, allowing for a more complete nuclear network evolution.
- **Python Bindings**: We plan to add Python bindings to GridFire, allowing for easy integration with Python-based workflows.

### Low Priority
- **GPU and Multi-Core Support**: We plan to add support for GPU and multi-core execution, allowing for faster network evolution.
- **Topology Optimization**: We plan to add support for topology optimization, allowing for more efficient network evolution.
- **Real-Time Network Visualization**: We plan to add support for real-time network visualization, allowing for better understanding of the network evolution.

## Building
GridFire uses meson as its build system. The minimum required version is 1.5.0. To build GridFire, you will need to have the following dependencies installed:

- C++ compiler supporting at least C++20 (though we test against C++23, and it is **strongly** recommended to use C++23)
- [Meson](https://mesonbuild.com/) build system (`apt install meson` or `brew install meson` or `pip install "meson>=1.5.0"`)
- [Ninja](https://ninja-build.org/) build system (`apt install ninja-build` or `brew install ninja` or `pip install ninja`)
- (_optional to have preinstalled_) [Boost](https://www.boost.org/) libraries (`apt install libboost-all-dev` or `brew install boost`)

Boost is labeled as optional because if GridFire cannot find boost during installation it will ask the user if they want
to install it, and then it will take care of that insallation. 

Once at least the C++ compiler, meson, and ninja are installed, you can build GridFire by running the following commands in the root directory of the repository:

```bash
meson setup build --buildtype=release
meson compile -C build
meson install -C build
```

Running the first command will create a `build` directory in the root of the repository, which will contain all the
build files. The second command will compile the library and the third command will install it to your system.
When installing GridFire it will generate a `pkg-config` file which can be used to link against the library in other
projects.

### Linking Against GridFire
The easiest way, by far, to link against GridFire is to use `pkg-config`. This will automatically handle all the
library and include flags for you. Further, the `pkg-config` file will also include the required information on
the libcomposition library so that you do not need to also manually link against it.

```bash
g++ -o main.cpp $(pkg-config --cflags --libs gridfire)
```

### Meson Wrap Installed Dependencies
Aside from Boost, GridFire dependes on a few other libraries which are not required to be installed on the system, but
are instead downloaded and built as part of the GridFire build process. These are:
- [CppAD](https://github.com/coin-or/CppAD): Used for automatic differentiation.
- [libcomposition](https://github.com/4d-star/libcomposition): Used for composition of networks.
- [libconfig](https://github.com/4d-star/libconfig): Used for configuration management.
- [libconstants](https://github.com/4d-star/libconstants): Used for physical constants.
- [liblogging](https://github.com/4d-star/liblogging): Used for logging.

These will all be downloaded and build automatically by meson when you run the `meson setup` and `meson compile` commands.

## Usage

> All code in GridFire is within the `gridfire` namespace. This means that you will need to use `gridfire::` prefix when
> using any of the classes or functions in GridFire.

> All code in lib* libraries is within the `fourdst::` namespace or subnamespaces of this (`fourdst::composition`,
`fourdst::config`, etc...). This means that you will need to use `fourdst::` prefix when using any of the classes or
> functions in these libraries.

GridFire is designed to primarily interface with `fourdst::composition::Composition` objects (provided by libcomposition).
The idea is that you create a composition object which tracks some set of species and either their mass fractions or number
fractions. You then initialize a network from this composition object. 

When evaluating the network you pass the current state of the composition object, the temperature, the density, and the 
current energy of the system (along with some other parmeters such as the maximum time and the initial time step). Any
network implemented in GridFire will then use this information to evolve the network and return a `gridfire::NetOut` object
which contains the updated state of the composition object, the new energy of the system, and the number of steps the
solver took to reach the final state.

There are a few networks currently implemented in GridFire.

| Network Name               | Description                                                                                                                                                                             | Network Format Enum |
|----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------|
| `gridfire::GraphNetwork`   | A general nuclear network which uses a directed hypergraph to represent the network topology. This is the main network implemented in GridFire and is the one you will most likely use. | `REACLIB`         |
| `gridfire::StaticGraphNetwork` |  A less general version of `gridfire::GraphNetwork` which uses a static graph to represent the network topology. That is to say that new species cannot be added or removed after instatiation of the network. | `REACLIB_STATIC`  |
| `gridfire::Approx8Network` | A specialized network for the approximate 8 species network used in stellar nucleosynthesis. This is a reimplimentation of the network developed by Frank Timmes' for MESA              | `APPROX8`          |
| `gridfire::Network`        | A virtual base class for all networks in GridFire. This is not meant to be used directly; rather, all networks in GridFire should subclass `gridfire::Network` and implement the required `evaluate` method. | `DEFAULT`         |

Regardless of the network you chose, the interface is the same. You will need to create a `gridfire::NetIn` object
which contains the initial state of the network, including the composition, temperature, density, energy, initial time
step, and maximum time. You will then pass this object to the `evaluate` method of the network you chose. The `evaluate`
method will return a `gridfire::NetOut` object which contains the updated state of the network.

A simple example of how to use `gridfire::GraphNetwork` is shown below (note that `composition.h` is a header file
provided by libcomposition, which means your compiler must be able to find its header and library files. The pkg-config
file generated during the installation of GridFire will take care of this for you, so you can use `pkg-config` to
compile this example):

```c++
#include <gridfire/netgraph.h>
#include <fourdst/composition/composition.h>

#include <iostream>
#include <vector>
#include <string>

int main() {
    const std::vector<double> comp = {0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    const std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};

    fourdst::composition::Composition composition;
    
    // true puts this in mass fraction mode. false would mean you must provide number fractions instead.
    composition.registerSymbol(symbols, true); 
    
    composition.setMassFraction(symbols, comp);
    composition.finalize(true);

    gridfire::NetIn netIn;
    netIn.composition = composition;
    netIn.temperature = 1e7; // In Kelvin
    netIn.density = 1e2;     // In g/cm^3
    netIn.energy = 0.0;      // Initial energy in erg

    netIn.tMax = 3.15e17;    // In seconds
    netIn.dt0 = 1e12;        // Initial time step in seconds

    gridfire::GraphNetwork network(composition);
    network.exportToDot("Network.dot");
    gridfire::NetOut netOut;
    netOut = network.evaluate(netIn);
    std::cout << netOut << std::endl;    
}

```
Save that file to `main.cpp` and compile it with the following command:
```bash
g++ main.cpp $(pkg-config --cflags --libs gridfire) -o main
./main
```

Using a different network is as simple as changing the type of the `network` variable to the desired network type. For
example, if you wanted to use the `gridfire::Approx8Network`, you would change the line
`gridfire::GraphNetwork network(composition);` to `gridfire::Approx8Network network(composition);`. The rest of the code
would remain the same.

For more details on how to use the composition library please refer to the libcomposition documentation.

## Funding
GridFire is a part of the 4D-STAR collaboration.

4D-STAR is funded by European Research Council (ERC) under the Horizon Europe programme (Synergy Grant agreement No.
101071505: 4D-STAR)
Work for this project is funded by the European Union. Views and opinions expressed are however those of the author(s)
only and do not necessarily reflect those of the European Union or the European Research.