#pragma once

#include "gridfire/network.h" // For NetIn, NetOut
#include "../reaction/reaction.h"

#include "fourdst/composition/composition.h"
#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"

#include <vector>
#include <unordered_map>

namespace gridfire {

    template<typename T>
    concept IsArithmeticOrAD = std::is_same_v<T, double> || std::is_same_v<T, CppAD::AD<double>>;

    template <IsArithmeticOrAD T>
    struct StepDerivatives {
        std::vector<T> dydt; ///< Derivatives of abundances.
        T nuclearEnergyGenerationRate = T(0.0); ///< Specific energy generation rate.
    };

    class Engine {
    public:
        virtual ~Engine() = default;
        virtual const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const = 0;
        virtual StepDerivatives<double> calculateRHSAndEnergy(
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const = 0;
    };

    class DynamicEngine : public Engine {
    public:
        virtual void generateJacobianMatrix(
            const std::vector<double>& Y,
            double T9, double rho
        ) = 0;
        virtual double getJacobianMatrixEntry(
            int i,
            int j
        ) const = 0;

        virtual void generateStoichiometryMatrix() = 0;
        virtual int getStoichiometryMatrixEntry(
            int speciesIndex,
            int reactionIndex
        ) const = 0;

        virtual double calculateMolarReactionFlow(
            const reaction::Reaction& reaction,
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const = 0;

        virtual const reaction::REACLIBLogicalReactionSet& getNetworkReactions() const = 0;
        virtual std::unordered_map<fourdst::atomic::Species, double> getSpeciesTimescales(
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const = 0;
    };
}