#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/exceptions/error_solver.h"
#include "gridfire/types/types.h"

#include <functional>
#include <any>
#include <vector>
#include <tuple>
#include <string>

namespace gridfire::solver {
    template <typename EngineT>
    concept IsEngine = std::is_base_of_v<engine::Engine, EngineT>;

    struct SolverContextBase {
        virtual void init() = 0;
        virtual void set_stdout_logging(bool enable) = 0;
        virtual void set_detailed_logging(bool enable) = 0;
        virtual ~SolverContextBase() = default;
    };

    /**
     * @struct TimestepContextBase
     * @brief Base class for solver callback contexts.
     *
     * This struct serves as a base class for contexts that can be passed to solver callbacks, it enforces
     * that derived classes implement a `describe` method that returns a vector of tuples describing
     * the context that a callback will receive when called.
     */
    class TimestepContextBase {
    public:
        virtual ~TimestepContextBase() = default;

        /**
         * @brief Describe the context for callback functions.
         * @return A vector of tuples, each containing a string for the parameters name and a string for its type.
         *
         * This method should be overridden by derived classes to provide a description of the context
         * that will be passed to the callback function. The intent of this method is that an end user can investigate
         * the context that will be passed to the callback function, and use this information to craft their own
         * callback function.
         */
        [[nodiscard]] virtual std::vector<std::tuple<std::string, std::string>> describe() const = 0;
    };
    /**
     * @brief Instantaneous rates at netIn's state straight from the engine, with no integration.
     *
     * One calculateRHSAndEnergy call (eps, neutrino losses, dY/dt) and one calculateEpsDerivatives call, plus
     * calculateEnergyChannels when include_channels is set. netIn.tMax, dt0 and energy are ignored. Network species
     * missing from netIn.composition are taken as absent.
     *
     * @param engine_ctx Scratch state of engine (its own, or one cloned from it); only engine scratch is touched.
     * @throws exceptions::SolverError for a non-positive or non-finite temperature or density, a negative or
     *         non-finite abundance, or a positive abundance of a species that is not in the network.
     * @throws exceptions::BadRHSEngineError if the engine cannot evaluate the RHS.
     */
    [[nodiscard]] InstantaneousRates evaluate_instantaneous_rates(
        const engine::DynamicEngine& engine,
        engine::scratch::StateBlob& engine_ctx,
        const NetIn& netIn,
        bool include_channels = true
    );

    /**
     * @class SingleZoneNetworkSolver
     * @brief Abstract base class for network solver strategies.
     *
     * This class defines the interface for network solver strategies, which are responsible
     * for integrating the reaction network ODEs over a given timestep.  It is templated on the
     * engine type to allow for different engine implementations to be used with the same solver.
     *
     * @tparam EngineT The type of engine to use with this solver strategy.  Must inherit from Engine.
     */
    template <IsEngine EngineT>
    class SingleZoneNetworkSolver {
    public:
        /**
         * @brief Constructor for the NetworkSolverStrategy.
         * @param engine The engine to use for evaluating the network.
         */
        explicit SingleZoneNetworkSolver(
            const EngineT& engine
        ) :
        m_engine(engine) {};

        /**
         * @brief Virtual destructor.
         */
        virtual ~SingleZoneNetworkSolver() = default;

        /**
         * @brief Evaluates the network for a given timestep.
         * @param solver_ctx
         * @param engine_ctx
         * @param netIn The input conditions for the network.
         * @return The output conditions after the timestep.
         */
        virtual NetOut evaluate(
            SolverContextBase& solver_ctx,
            const NetIn& netIn
        ) const = 0;

        /**
         * @brief Instantaneous eps, its T and rho derivatives, neutrino losses, dY/dt and channel energies at netIn's
         *        state, without integrating (netIn.tMax and dt0 are ignored).
         *
         * A pass-through to the engine (see evaluate_instantaneous_rates) using this solver's context, so a consumer
         * holding a solver stack can get the rates that evaluate() would report at the end of a step, at any state.
         * The integrator state in solver_ctx is not touched. Solvers that support it override this; the default
         * throws exceptions::SolverError.
         */
        [[nodiscard]] virtual InstantaneousRates evaluateRates(
            SolverContextBase& solver_ctx,
            const NetIn& netIn
        ) const {
            throw exceptions::SolverError("This solver does not support evaluateRates.");
        }

        /**
         * @brief Burn one composition shared by several shells, each at its own temperature and density, for
         *        netIn.tMax: dY/dt = sum_i w_i f(Y, T_i, rho_i). netIn supplies the composition, tMax and dt0; its
         *        temperature and density are ignored.
         *
         * This is a convective zone under instantaneous mixing: the shells exchange material much faster than they
         * burn, so they hold one composition and the zone burns at the mass-weighted mean of its shells' rates. Each
         * shell's energy and neutrino losses are integrated separately. The default throws exceptions::SolverError.
         *
         * @param shells The shells' states and mass weights (weights summing to 1).
         */
        [[nodiscard]] virtual SharedCompositionResult evaluateSharedComposition(
            SolverContextBase& solver_ctx,
            const NetIn& netIn,
            const std::vector<ShellState>& shells
        ) const {
            throw exceptions::SolverError("This solver does not support evaluateSharedComposition.");
        }

    protected:
        const EngineT& m_engine; ///< The engine used by this solver strategy.
    };

    template <IsEngine EngineT>
    class MultiZoneNetworkSolver {
    public:
        explicit MultiZoneNetworkSolver(
            const EngineT& engine,
            const SingleZoneNetworkSolver<EngineT>& solver
        ) :
        m_engine(engine),
        m_solver(solver) {};

        virtual ~MultiZoneNetworkSolver() = default;

        virtual std::vector<NetOut> evaluate(
            SolverContextBase& solver_ctx,
            const std::vector<NetIn>& netIns
        ) const = 0;

        /**
         * @brief SingleZoneNetworkSolver::evaluateRates for every zone, without integrating. The default throws
         *        exceptions::SolverError.
         */
        [[nodiscard]] virtual std::vector<InstantaneousRates> evaluateRates(
            SolverContextBase& solver_ctx,
            const std::vector<NetIn>& netIns
        ) const {
            throw exceptions::SolverError("This solver does not support evaluateRates.");
        }
    protected:
        const EngineT& m_engine; ///< The engine used by this solver strategy.
        const SingleZoneNetworkSolver<EngineT>& m_solver;
    };

    /**
     * @brief Type alias for a network solver strategy that uses a DynamicEngine.
     */
    using SingleZoneDynamicNetworkSolver = SingleZoneNetworkSolver<engine::DynamicEngine>;
    using MultiZoneDynamicNetworkSolver = MultiZoneNetworkSolver<engine::DynamicEngine>;
}
