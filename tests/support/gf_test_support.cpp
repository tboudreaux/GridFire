#include "gf_test_support.h"

#include "fourdst/composition/utils/utils.h"
#include "fourdst/constants/const.h"
#include "fourdst/logging/logging.h"
#include "gridfire/utils/gf_omp.h"

#include "quill/Logger.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <format>
#include <limits>
#include <stdexcept>

namespace gridfire::test {

    void quiet_logging() {
        fourdst::logging::LogManager::getInstance().getLogger("log")->set_log_level(quill::LogLevel::Warning);
    }

    fourdst::composition::Composition solar_composition() {
        const std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
        const std::vector<double> X = {0.7081145999999999, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
        return fourdst::composition::buildCompositionFromMassFractions(symbols, X);
    }

    namespace {
        struct EngineHolder {
            std::unique_ptr<policy::MainSequencePolicy> policy;
            std::unique_ptr<policy::ConstructionResults> results;
            std::unique_ptr<SharedEngine> shared;
        };

        EngineHolder& holder() {
            static EngineHolder h = [] {
                quiet_logging();
                GF_PAR_INIT(); // Must precede any CppAD use so its allocator can tell OpenMP threads apart
                EngineHolder built;
                built.policy = std::make_unique<policy::MainSequencePolicy>(solar_composition());
                built.results = std::make_unique<policy::ConstructionResults>(built.policy->construct());
                const auto* graph = dynamic_cast<const engine::GraphEngine*>(&built.results->engine);
                if (!graph) {
                    throw std::runtime_error("MainSequencePolicy no longer returns a bare GraphEngine; update gf_test_support.");
                }
                built.shared = std::make_unique<SharedEngine>(SharedEngine{built.results->engine, *graph});
                return built;
            }();
            return h;
        }
    }

    const SharedEngine& main_sequence_engine() {
        return *holder().shared;
    }

    std::unique_ptr<engine::scratch::StateBlob> SharedEngine::fresh_context() const {
        return holder().results->scratch_blob->clone_structure();
    }

    std::string EngineState::describe() const {
        return std::format("T9={:.6e}, rho={:.6e}, species={}", T9, rho, composition.size());
    }

    std::uint64_t test_seed() {
        if (const char* env = std::getenv("GF_TEST_SEED")) {
            return std::strtoull(env, nullptr, 10);
        }
        return 20261104ULL;
    }

    StateSampler::StateSampler(
        const std::vector<fourdst::atomic::Species>& species,
        const std::uint64_t seed,
        const double T9_min,
        const double T9_max,
        const double rho_min,
        const double rho_max
    ) : m_species(species), m_rng(seed), m_T9_min(T9_min), m_T9_max(T9_max), m_rho_min(rho_min), m_rho_max(rho_max) {}

    EngineState StateSampler::next() {
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        auto log_uniform = [&](const double lo, const double hi) {
            return std::exp(std::log(lo) + unit(m_rng) * (std::log(hi) - std::log(lo)));
        };

        const auto solar = solar_composition();
        std::vector<double> X(m_species.size());
        for (size_t i = 0; i < m_species.size(); ++i) {
            const auto& sp = m_species[i];
            if (solar.contains(sp)) {
                X[i] = solar.getMassFraction(sp) * log_uniform(0.1, 10.0);
            } else {
                X[i] = log_uniform(1e-20, 1e-3);
            }
        }
        double total = 0.0;
        for (const double x : X) total += x;

        std::vector<double> Y(m_species.size());
        for (size_t i = 0; i < m_species.size(); ++i) {
            Y[i] = (X[i] / total) / m_species[i].mass();
        }

        return {log_uniform(m_T9_min, m_T9_max), log_uniform(m_rho_min, m_rho_max), make_composition(m_species, Y)};
    }

    fourdst::composition::Composition make_composition(
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<double>& Y
    ) {
        return {species, Y};
    }

    std::vector<double> abundances_in_order(
        const fourdst::composition::CompositionAbstract& comp,
        const std::vector<fourdst::atomic::Species>& species
    ) {
        std::vector<double> Y;
        Y.reserve(species.size());
        for (const auto& sp : species) {
            Y.push_back(comp.contains(sp) ? comp.getMolarAbundance(sp) : 0.0);
        }
        return Y;
    }

    double rel_diff(const double a, const double b, const double floor) {
        const double scale = std::max({std::abs(a), std::abs(b), floor});
        if (scale == 0.0) return 0.0;
        return std::abs(a - b) / scale;
    }

    FiniteDifference richardson_derivative(const std::function<double(double)>& f, const double x, const double h) {
        const double fp = f(x + h), fm = f(x - h);
        const double fph = f(x + h / 2), fmh = f(x - h / 2);
        const double d_h = (fp - fm) / (2 * h);
        const double d_h2 = (fph - fmh) / h;
        const double fmax = std::max({std::abs(fp), std::abs(fm), std::abs(fph), std::abs(fmh)});
        return {(4 * d_h2 - d_h) / 3, 64 * std::numeric_limits<double>::epsilon() * fmax / h};
    }

    const PhysicalConstants& constants() {
        static const PhysicalConstants c = [] {
            const auto& k = fourdst::constant::Constants::getInstance();
            return PhysicalConstants{k.get("u").value, k.get("N_a").value, k.get("c").value};
        }();
        return c;
    }
}
