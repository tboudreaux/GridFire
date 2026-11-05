#include "gridfire/solver/strategies/strategy_abstract.h"

#include "gridfire/exceptions/error_engine.h"
#include "gridfire/exceptions/error_solver.h"

#include "fourdst/composition/composition.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <vector>

namespace gridfire::solver {
    InstantaneousRates evaluate_instantaneous_rates(
        const engine::DynamicEngine& engine,
        engine::scratch::StateBlob& engine_ctx,
        const NetIn& netIn,
        const bool include_channels
    ) {
        if (!std::isfinite(netIn.temperature) || netIn.temperature <= 0.0) {
            throw exceptions::SolverError(std::format("evaluateRates: temperature must be finite and positive (got {} K)", netIn.temperature));
        }
        if (!std::isfinite(netIn.density) || netIn.density <= 0.0) {
            throw exceptions::SolverError(std::format("evaluateRates: density must be finite and positive (got {} g/cm^3)", netIn.density));
        }

        // The engine needs every network species; those missing from the input are absent
        const auto& network = engine.getNetworkSpecies(engine_ctx);
        std::vector<double> Y(network.size(), 0.0);
        size_t matched = 0;
        for (size_t k = 0; k < network.size(); ++k) {
            if (netIn.composition.contains(network[k])) {
                Y[k] = netIn.composition.getMolarAbundance(network[k]);
                ++matched;
            }
            if (!std::isfinite(Y[k]) || Y[k] < 0.0) {
                throw exceptions::SolverError(std::format(
                    "evaluateRates: molar abundance of {} must be finite and non-negative (got {})", network[k].name(), Y[k]));
            }
        }
        if (matched != netIn.composition.size()) {
            for (const auto& species : netIn.composition.getRegisteredSpecies()) {
                if (std::ranges::find(network, species) == network.end() && netIn.composition.getMolarAbundance(species) > 0.0) {
                    throw exceptions::SolverError(std::format(
                        "evaluateRates: {} (molar abundance {:.3e}) is not in the network.", species.name(),
                        netIn.composition.getMolarAbundance(species)));
                }
            }
        }
        const fourdst::composition::Composition comp(network, Y);

        const double T9 = netIn.temperature / 1e9;
        const auto rhs = engine.calculateRHSAndEnergy(engine_ctx, comp, T9, netIn.density, false);
        if (!rhs) {
            throw exceptions::BadRHSEngineError(std::format(
                "evaluateRates: the engine could not evaluate the RHS ({}).", engine::EngineStatus_to_string(rhs.error())));
        }
        const auto [dEps_dT, dEps_dRho] = engine.calculateEpsDerivatives(engine_ctx, comp, T9, netIn.density);

        InstantaneousRates rates;
        rates.eps_nuc = rhs->nuclearEnergyGenerationRate;
        rates.dEps_dT = dEps_dT;
        rates.dEps_dRho = dEps_dRho;
        rates.neutrino_energy_loss_rate = rhs->neutrinoEnergyLossRate;
        rates.neutrino_flux = rhs->totalNeutrinoFlux;
        rates.dydt = rhs->dydt;
        if (include_channels) {
            rates.energy_channel_names = engine.getEnergyChannelNames(engine_ctx);
            rates.eps_by_channel = engine.calculateEnergyChannels(engine_ctx, comp, T9, netIn.density);
        }
        return rates;
    }
}
