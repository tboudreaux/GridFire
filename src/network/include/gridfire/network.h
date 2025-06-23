/* ***********************************************************************
//
//   Copyright (C) 2025 -- The 4D-STAR Collaboration
//   File Authors: Emily Boudreaux, Aaron Dotter
//   Last Modified: March 21, 2025
//
//   4DSSE is free software; you can use it and/or modify
//   it under the terms and restrictions the GNU General Library Public
//   License version 3 (GPLv3) as published by the Free Software Foundation.
//
//   4DSSE is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Library General Public License for more details.
//
//   You should have received a copy of the GNU Library General Public License
//   along with this software; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// *********************************************************************** */
#pragma once

#include <vector>

#include "fourdst/logging/logging.h"
#include "fourdst/config/config.h"
#include "fourdst/composition/species.h"
#include "quill/Logger.h"
#include "fourdst/composition/composition.h"
#include "gridfire/reaclib.h"
#include <unordered_map>

#include "fourdst/constants/const.h"

namespace gridfire {


    enum NetworkFormat {
        APPROX8, ///< Approx8 nuclear reaction network format.
        REACLIB, ///< General REACLIB nuclear reaction network format.
        UNKNOWN,
    };

    static inline std::unordered_map<NetworkFormat, std::string> FormatStringLookup = {
        {APPROX8, "Approx8"},
        {REACLIB, "REACLIB"},
        {UNKNOWN, "Unknown"}
    };

    /**
     * @struct NetIn
     * @brief Input structure for the network evaluation.
     *
     * This structure holds the input parameters required for the network evaluation.
     *
     * Example usage:
     * @code
     * nuclearNetwork::NetIn netIn;
     * netIn.composition = {1.0, 0.0, 0.0};
     * netIn.tmax = 1.0e6;
     * netIn.dt0 = 1.0e-3;
     * netIn.temperature = 1.0e8;
     * netIn.density = 1.0e5;
     * netIn.energy = 1.0e12;
     * @endcode
     */
    struct NetIn {
        fourdst::composition::Composition composition; ///< Composition of the network
        double tMax; ///< Maximum time
        double dt0; ///< Initial time step
        double temperature; ///< Temperature in Kelvin
        double density; ///< Density in g/cm^3
        double energy; ///< Energy in ergs
        double culling = 0.0; ///< Culling threshold for reactions (default is 0.0, meaning no culling)
    };

    /**
     * @struct NetOut
     * @brief Output structure for the network evaluation.
     * 
     * This structure holds the output results from the network evaluation.
     * 
     * Example usage:
     * @code
     * nuclearNetwork::NetOut netOut = network.evaluate(netIn);
     * std::vector<double> composition = netOut.composition;
     * int steps = netOut.num_steps;
     * double energy = netOut.energy;
     * @endcode
     */
    struct NetOut {
        fourdst::composition::Composition composition; ///< Composition of the network after evaluation
        int num_steps; ///< Number of steps taken in the evaluation
        double energy; ///< Energy in ergs after evaluation

        friend std::ostream& operator<<(std::ostream& os, const NetOut& netOut) {
            os << "NetOut(composition=" << netOut.composition << ", num_steps=" << netOut.num_steps << ", energy=" << netOut.energy << ")";
            return os;
        }
    };

    /**
     * @class Network
     * @brief Class for network evaluation.
     * 
     * This class provides methods to evaluate the network based on the input parameters.
     * 
     * Example usage:
     * @code
     * nuclearNetwork::Network network;
     * nuclearNetwork::NetIn netIn;
     * // Set netIn parameters...
     * nuclearNetwork::NetOut netOut = network.evaluate(netIn);
     * @endcode
     */
    class Network {
        public:
            explicit Network(const NetworkFormat format = NetworkFormat::APPROX8);
            virtual ~Network() = default;

            [[nodiscard]] NetworkFormat getFormat() const;
            NetworkFormat setFormat(const NetworkFormat format);

            /**
             * @brief Evaluate the network based on the input parameters.
             * 
             * @param netIn Input parameters for the network evaluation.
             * @return NetOut Output results from the network evaluation.
             */
            virtual NetOut evaluate(const NetIn &netIn);

            virtual bool isStiff() const { return m_stiff; }
            virtual void setStiff(const bool stiff) { m_stiff = stiff; }

        protected:
            fourdst::config::Config& m_config; ///< Configuration instance
            fourdst::logging::LogManager& m_logManager; ///< Log manager instance
            quill::Logger* m_logger; ///< Logger instance

            NetworkFormat m_format; ///< Format of the network
            fourdst::constant::Constants& m_constants;

            bool m_stiff = false; ///< Flag indicating if the network is stiff
    };

    class LogicalReaction {
    public:
        explicit LogicalReaction(const std::vector<reaclib::REACLIBReaction> &reactions);
        explicit LogicalReaction(const reaclib::REACLIBReaction &reaction);
        void add_reaction(const reaclib::REACLIBReaction& reaction);

        template <typename T>
        [[nodiscard]] T calculate_rate(const T T9) const {
            T sum = static_cast<T>(0.0);
            const T T913 = CppAD::pow(T9, 1.0/3.0);
            const T T953 = CppAD::pow(T9, 5.0/3.0);
            const T logT9 = CppAD::log(T9);
            for (const auto& rate : m_rates) {
                const T exponent = rate.a0 +
                       rate.a1 / T9 +
                       rate.a2 / T913 +
                       rate.a3 * T913 +
                       rate.a4 * T9 +
                       rate.a5 * T953 +
                       rate.a6 * logT9;
                sum += CppAD::exp(exponent);
            }
            return sum;
        }

        [[nodiscard]] const std::string& id() const {return std::string(m_peID); }

        [[nodiscard]] std::string_view peName() const { return m_peID; }

        [[nodiscard]] int chapter() const { return m_chapter; }

        [[nodiscard]] const std::vector<fourdst::atomic::Species>& reactants() const { return m_reactants; }

        [[nodiscard]] const std::vector<fourdst::atomic::Species>& products() const { return m_products; }

        [[nodiscard]] double qValue() const { return m_qValue; }

        [[nodiscard]] bool is_reverse() const { return m_reverse; }


        auto begin();
        auto begin() const;
        auto end();
        auto end() const;

    private:
        const quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::string_view m_peID;
        std::vector<fourdst::atomic::Species> m_reactants; ///< Reactants of the reaction
        std::vector<fourdst::atomic::Species> m_products; ///< Products of the reaction
        double m_qValue = 0.0; ///< Q-value of the reaction
        bool m_reverse = false; ///< True if the reaction is reverse
        int m_chapter;

        std::vector<reaclib::RateFitSet> m_rates;

    };

    class LogicalReactionSet {
    public:
        LogicalReactionSet() = default;
        explicit LogicalReactionSet(const std::vector<LogicalReaction>& reactions);
        explicit LogicalReactionSet(const std::vector<reaclib::REACLIBReaction>& reactions);
        explicit LogicalReactionSet(const reaclib::REACLIBReactionSet& reactionSet);

        void add_reaction(const LogicalReaction& reaction);
        void add_reaction(const reaclib::REACLIBReaction& reaction);

        void remove_reaction(const LogicalReaction& reaction);

        [[nodiscard]] bool contains(const std::string_view& id) const;
        [[nodiscard]] bool contains(const LogicalReaction& reactions) const;
        [[nodiscard]] bool contains(const reaclib::REACLIBReaction& reaction) const;

        [[nodiscard]] size_t size() const;

        void sort(double T9=1.0);

        bool contains_species(const fourdst::atomic::Species &species) const;
        bool contains_reactant(const fourdst::atomic::Species &species) const;
        bool contains_product(const fourdst::atomic::Species &species) const;

        [[nodiscard]] const LogicalReaction& operator[](size_t index) const;
        [[nodiscard]] const LogicalReaction& operator[](const std::string_view& id) const;

        auto begin();
        auto begin() const;
        auto end();
        auto end() const;


    private:
        const quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::vector<LogicalReaction> m_reactions;
        std::unordered_map<std::string_view, LogicalReaction&> m_reactionNameMap;
    };

    LogicalReactionSet build_reaclib_nuclear_network(const fourdst::composition::Composition &composition);
    LogicalReactionSet build_reaclib_nuclear_network(const fourdst::composition::Composition &composition, double culling, double T9 = 1.0);


} // namespace nuclearNetwork
