#include <gtest/gtest.h>
#include <string>

#include "approx8.h"
#include "config.h"
#include "network.h"

#include <vector>

std::string TEST_CONFIG = std::string(getenv("MESON_SOURCE_ROOT")) + "/tests/testsConfig.yaml";
class approx8Test : public ::testing::Test {};

/**
 * @brief Test the constructor of the Config class.
 */
TEST_F(approx8Test, constructor) {
    Config& config = Config::getInstance();
    config.loadConfig(TEST_CONFIG);
    EXPECT_NO_THROW(nnApprox8::Approx8Network());
}

TEST_F(approx8Test, setStiff) {
    nnApprox8::Approx8Network network;
    EXPECT_NO_THROW(network.setStiff(true));
    EXPECT_TRUE(network.isStiff());
    EXPECT_NO_THROW(network.setStiff(false));
    EXPECT_FALSE(network.isStiff());
}

TEST_F(approx8Test, evaluate) {
    nnApprox8::Approx8Network network;
    nuclearNetwork::NetIn netIn;
    
    std::vector<double> comp = {0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};

    netIn.composition = comp;
    netIn.temperature = 1e7;
    netIn.density = 1e2;
    netIn.energy = 0.0;
    
    netIn.tmax = 3.15e17;
    netIn.dt0 = 1e12;

    nuclearNetwork::NetOut netOut;
    EXPECT_NO_THROW(netOut = network.evaluate(netIn));

    EXPECT_DOUBLE_EQ(netOut.composition[nnApprox8::Net::ih1], 0.50166260916650918);
    EXPECT_DOUBLE_EQ(netOut.composition[nnApprox8::Net::ihe4],0.48172270591286032);
    EXPECT_DOUBLE_EQ(netOut.energy, 1.6433049870528356e+18);
}
