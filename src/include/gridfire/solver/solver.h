#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"

#include <functional>
#include <any>
#include <vector>
#include <tuple>
#include <string>

namespace gridfire::solver {
    /**
     * @struct SolverContextBase
     * @brief Base class for solver callback contexts.
     *
     * This struct serves as a base class for contexts that can be passed to solver callbacks, it enforces
     * that derived classes implement a `describe` method that returns a vector of tuples describing
     * the context that a callback will receive when called.
     */
    struct SolverContextBase {
        virtual ~SolverContextBase() = default;

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
     * @class NetworkSolverStrategy
     * @brief Abstract base class for network solver strategies.
     *
     * This class defines the interface for network solver strategies, which are responsible
     * for integrating the reaction network ODEs over a given timestep.  It is templated on the
     * engine type to allow for different engine implementations to be used with the same solver.
     *
     * @tparam EngineT The type of engine to use with this solver strategy.  Must inherit from Engine.
     */
    template <typename EngineT>
    class NetworkSolverStrategy {
    public:
        /**
         * @brief Constructor for the NetworkSolverStrategy.
         * @param engine The engine to use for evaluating the network.
         */
        explicit NetworkSolverStrategy(EngineT& engine) : m_engine(engine) {};

        /**
         * @brief Virtual destructor.
         */
        virtual ~NetworkSolverStrategy() = default;

        /**
         * @brief Evaluates the network for a given timestep.
         * @param netIn The input conditions for the network.
         * @return The output conditions after the timestep.
         */
        virtual NetOut evaluate(const NetIn& netIn) = 0;

        /**
         * @brief set the callback function to be called at the end of each timestep.
         *
         * This function allows the user to set a callback function that will be called at the end of each timestep.
         * The callback function will receive a gridfire::solver::<SOMESOLVER>::TimestepContext object. Note that
         * depending on the solver, this context may contain different information. Further, the exact
         * signature of the callback function is left up to each solver. Every solver should provide a type or type alias
         * TimestepCallback that defines the signature of the callback function so that the user can easily
         * get that type information.
         *
         * @param callback The callback function to be called at the end of each timestep.
         */
        virtual void set_callback(const std::any& callback) = 0;

        /**
         * @brief Describe the context that will be passed to the callback function.
         * @return A vector of tuples, each containing a string for the parameter's name and a string for its type.
         *
         * This method should be overridden by derived classes to provide a description of the context
         * that will be passed to the callback function. The intent of this method is that an end user can investigate
         * the context that will be passed to the callback function, and use this information to craft their own
         * callback function.
         */
        [[nodiscard]] virtual std::vector<std::tuple<std::string, std::string>> describe_callback_context() const = 0;
    protected:
        EngineT& m_engine; ///< The engine used by this solver strategy.
    };

    /**
     * @brief Type alias for a network solver strategy that uses a DynamicEngine.
     */
    using DynamicNetworkSolverStrategy = NetworkSolverStrategy<DynamicEngine>;
}