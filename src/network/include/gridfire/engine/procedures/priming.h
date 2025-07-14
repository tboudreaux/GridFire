#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"

#include "fourdst/composition/composition.h"
#include "fourdst/composition/atomicSpecies.h"

#include <map>
#include <ranges>
#include <sstream>


namespace gridfire {


    PrimingReport primeNetwork(
        const NetIn&,
        DynamicEngine& engine
    );

    double calculateDestructionRateConstant(
        const DynamicEngine& engine,
        const fourdst::atomic::Species& species,
        const std::vector<double>& Y,
        double T9,
        double rho
    );

    double calculateCreationRate(
        const DynamicEngine& engine,
        const fourdst::atomic::Species& species,
        const std::vector<double>& Y,
        double T9,
        double rho
    );
}