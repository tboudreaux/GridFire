#include "gridfire/engine/procedures/priming.h"

#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils/utils.h"
#include "gridfire/solver/solver.h"

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/types/types.h"
#include "gridfire/exceptions/error_solver.h"

#include "gridfire/engine/scratchpads/blob.h"
#include "gridfire/engine/scratchpads/engine_graph_scratchpad.h"

#include "fourdst/logging/logging.h"
#include "gridfire/solver/strategies/PointSolver.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"


namespace gridfire::engine {
    using fourdst::composition::Composition;
    using fourdst::atomic::Species;

    PrimingReport primeNetwork(
        scratch::StateBlob &ctx,
        const NetIn& netIn,
        const GraphEngine& engine, const std::optional<std::vector<reaction::ReactionType>>& ignoredReactionTypes
    ) {
        const auto logger = LogManager::getInstance().getLogger("log");
        solver::PointSolver integrator(engine);
        solver::PointSolverContext solverCtx(ctx);
        solverCtx.abs_tol = 1e-3;
        solverCtx.rel_tol = 1e-3;
        solverCtx.stdout_logging = false;

        NetIn solverInput(netIn);

        solverInput.tMax = 1e-15;
        solverInput.temperature = 1e7;

        LOG_INFO(logger, "Short timescale ({}) network ignition started.", solverInput.tMax);
        PrimingReport report;
        try {
            const NetOut netOut = integrator.evaluate(solverCtx, solverInput);
            LOG_INFO(logger, "Network ignition completed.");
            LOG_TRACE_L2(
                logger,
                "After ignition composition is {}",
                [netOut, netIn]() -> std::string {
                    std::stringstream ss;
                    size_t i = 0;
                    for (const auto& [species, abundance] : netOut.composition) {
                        ss << species.name() << ": " << abundance << " (prior: " << netIn.composition.getMolarAbundance(species);
                        ss << ", fractional change: " << (abundance - netIn.composition.getMolarAbundance(species)) / netIn.composition.getMolarAbundance(species) * 100.0 << "%)";
                        if (i < netOut.composition.size() - 1) {
                            ss << ", ";
                        }
                        ++i;
                    }
                    return ss.str();
                }()
            );
            report.primedComposition = netOut.composition;
            std::unordered_set<Species> unprimedSpecies;
            double minAbundance = std::numeric_limits<double>::max();
            for (const auto& [sp, y] : report.primedComposition) {
                if (y == 0) {
                    unprimedSpecies.insert(sp);
                }
                if (y != 0 && y < minAbundance) {
                    minAbundance = y;
                }
            }
            const double abundanceForUnprimedSpecies = minAbundance / 1e10;
            for (const auto& sp : unprimedSpecies) {
                LOG_TRACE_L1(logger, "Clamping Species {}: initial abundance {}, primed abundance {} to {}", sp.name(), netIn.composition.getMolarAbundance(sp), report.primedComposition.getMolarAbundance(sp), abundanceForUnprimedSpecies);
                report.primedComposition.setMolarAbundance(sp, abundanceForUnprimedSpecies);
            }
            report.success = true;
            report.status = PrimingReportStatus::SUCCESS;
        } catch (const exceptions::SolverError& e) {
            LOG_ERROR(logger, "Failed to prime network: solver failure encountered: {}", e.what());
            std::rethrow_exception(std::current_exception());
        }
        return report;
    }
}