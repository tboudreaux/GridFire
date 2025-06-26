#define QUILL_COMPILE_ACTIVE_LOG_LEVEL QUILL_COMPILE_ACTIVE_LOG_LEVEL_DEBUG
#include <iostream>
#include "../../src/network/include/gridfire/engine/engine_graph.h"
#include "gridfire/network.h"
#include "../../src/network/include/gridfire/engine/engine_approx8.h"
#include "../../src/network/include/gridfire/solver/solver.h"

#include "fourdst/composition/composition.h"

int main() {
    using namespace gridfire;
    // reaction::REACLIBLogicalReactionSet reactions = build_reaclib_nuclear_network_from_file("approx8", false);
    // GraphNetwork network(reactions);
    const std::vector<double> comp = {0.708, 0.0, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    const std::vector<std::string> symbols = {"H-1", "H-2", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};

    fourdst::composition::Composition composition;
    composition.registerSymbol(symbols, true);
    composition.setMassFraction(symbols, comp);
    composition.finalize(true);


    NetIn netIn;
    netIn.composition = composition;
    netIn.temperature = 1e7;
    netIn.density = 1e2;
    netIn.energy = 0.0;

    netIn.tMax = 3.15e17;
    netIn.dt0 = 1e12;

    NetOut netOut;
    approx8::Approx8Network approx8Network;
    netOut = approx8Network.evaluate(netIn);
    std::cout << "Approx8 Network Output: " << netOut << std::endl;

    // // Network Ignition
    netIn.dt0 = 1e-15;
    // netIn.temperature = 2.0e8;
    // netIn.density = 1.0e6;
    // netIn.tMax = 1e-1;

    GraphEngine ReaclibEngine(composition);
    // netOut = network.evaluate(netIn);
    solver::QSENetworkSolver solver(ReaclibEngine);
    netOut = solver.evaluate(netIn);
    std::cout << "QSE Graph Network Output: " << netOut << std::endl;
}