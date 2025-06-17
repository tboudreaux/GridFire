#include <gtest/gtest.h>
#include <string>

#include "approx8.h"
#include "config.h"
#include "network.h"
#include "composition.h"

#include <vector>

std::string TEST_CONFIG = std::string(getenv("MESON_SOURCE_ROOT")) + "/tests/testsConfig.yaml";
class approx8Test : public ::testing::Test {};

/**
 * @brief Test the constructor of the Config class.
 */
TEST_F(approx8Test, constructor) {
    serif::config::Config& config = serif::config::Config::getInstance();
    config.loadConfig(TEST_CONFIG);
    EXPECT_NO_THROW(serif::network::approx8::Approx8Network());
}

TEST_F(approx8Test, setStiff) {
    serif::network::approx8::Approx8Network network;
    EXPECT_NO_THROW(network.setStiff(true));
    EXPECT_TRUE(network.isStiff());
    EXPECT_NO_THROW(network.setStiff(false));
    EXPECT_FALSE(network.isStiff());
}

TEST_F(approx8Test, evaluate) {
    serif::network::approx8::Approx8Network network;
    serif::network::NetIn netIn;

    std::vector<double> comp = {0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};

    serif::composition::Composition composition;
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

    serif::network::NetOut netOut;
    EXPECT_NO_THROW(netOut = network.evaluate(netIn));

    double energyFraction = netOut.energy / 1.6433051127589775E+18;
    double H1MassFraction = netOut.composition.getMassFraction("H-1")/ 0.50166262445895604;
    double He4MassFraction = netOut.composition.getMassFraction("He-4") / 0.48172273720971226;

    double relError = 1e-8;
    EXPECT_NEAR(H1MassFraction, 1.0, relError);
    EXPECT_NEAR(He4MassFraction, 1.0, relError);
    EXPECT_NEAR(energyFraction, 1.0, relError);
}
