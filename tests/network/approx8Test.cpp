#include <string>
#include <gtest/gtest.h>

#include "fourdst/composition/composition.h"
#include "fourdst/config/config.h"
#include "gridfire/engine/engine_approx8.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/network.h"

#include <vector>


std::string TEST_CONFIG = std::string(getenv("MESON_SOURCE_ROOT")) + "/tests/testsConfig.yaml";
class approx8Test : public ::testing::Test {};

/**
 * @brief Test the constructor of the Config class.
 */
TEST_F(approx8Test, constructor) {
    fourdst::config::Config& config = fourdst::config::Config::getInstance();
    config.loadConfig(TEST_CONFIG);
    EXPECT_NO_THROW(gridfire::approx8::Approx8Network());
}

TEST_F(approx8Test, setStiff) {
    gridfire::approx8::Approx8Network network;
    EXPECT_NO_THROW(network.setStiff(true));
    EXPECT_TRUE(network.isStiff());
    EXPECT_NO_THROW(network.setStiff(false));
    EXPECT_FALSE(network.isStiff());
}

TEST_F(approx8Test, evaluate) {
    gridfire::approx8::Approx8Network network;
    gridfire::NetIn netIn;

    std::vector<double> comp = {0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};

    fourdst::composition::Composition composition;
    composition.registerSymbol(symbols, true);
    composition.setMassFraction(symbols, comp);
    bool isFinalized = composition.finalize(true);
    EXPECT_TRUE(isFinalized);

    netIn.composition = composition;
    netIn.temperature = 1e7;
    netIn.density = 1e2;
    netIn.energy = 0.0;

    netIn.tMax = 3.15e17;
    netIn.dt0 = 1e12;

    gridfire::NetOut netOut;
    EXPECT_NO_THROW(netOut = network.evaluate(netIn));

    double energyFraction = netOut.energy / 1.6433051127589775E+18;
    double H1MassFraction = netOut.composition.getMassFraction("H-1")/ 0.50166262445895604;
    double He4MassFraction = netOut.composition.getMassFraction("He-4") / 0.48172273720971226;

    double relError = 1e-6;
    EXPECT_NEAR(H1MassFraction, 1.0, relError);
    EXPECT_NEAR(He4MassFraction, 1.0, relError);
    EXPECT_NEAR(energyFraction, 1.0, relError);
}

TEST_F(approx8Test, reaclib) {
    using namespace gridfire;
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

    GraphEngine network(composition);
    // NetOut netOut;
    // netOut = network.evaluate(netIn);
    // std::cout << netOut << std::endl;
}
