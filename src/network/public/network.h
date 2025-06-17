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

#include "probe.h"
#include "config.h"
#include "quill/Logger.h"
#include "composition.h"
#include <unordered_map>

namespace serif::network {

    enum NetworkFormat {
        APPROX8, ///< Approx8 nuclear reaction network format.
        UNKNOWN,
    };

    static inline std::unordered_map<NetworkFormat, std::string> FormatStringLookup = {
        {APPROX8, "Approx8"},
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
        serif::composition::Composition composition; ///< Composition of the network
        double tMax; ///< Maximum time
        double dt0; ///< Initial time step
        double temperature; ///< Temperature in Kelvin
        double density; ///< Density in g/cm^3
        double energy; ///< Energy in ergs
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
        serif::composition::Composition composition; ///< Composition of the network after evaluation
        int num_steps; ///< Number of steps taken in the evaluation
        double energy; ///< Energy in ergs after evaluation
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

            NetworkFormat getFormat() const;
            NetworkFormat setFormat(const NetworkFormat format);

            /**
             * @brief Evaluate the network based on the input parameters.
             * 
             * @param netIn Input parameters for the network evaluation.
             * @return NetOut Output results from the network evaluation.
             */
            virtual NetOut evaluate(const NetIn &netIn);

        protected:
            serif::config::Config& m_config; ///< Configuration instance
            serif::probe::LogManager& m_logManager; ///< Log manager instance
            quill::Logger* m_logger; ///< Logger instance

            NetworkFormat m_format; ///< Format of the network
    };



} // namespace nuclearNetwork
