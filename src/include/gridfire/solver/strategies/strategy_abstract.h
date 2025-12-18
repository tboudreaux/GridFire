#pragma once

#include "gridfire/engine/engine_abstract.h"
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
