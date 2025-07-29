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
#include "fourdst/composition/composition.h"
#include "fourdst/constants/const.h"

#include "gridfire/reaction/reaction.h"

#include "quill/Logger.h"

#include <unordered_map>


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

    struct NetIn {
        fourdst::composition::Composition composition; ///< Composition of the network
        double tMax; ///< Maximum time
        double dt0; ///< Initial time step
        double temperature; ///< Temperature in Kelvin
        double density; ///< Density in g/cm^3
        double energy; ///< Energy in ergs
        double culling = 0.0; ///< Culling threshold for reactions (default is 0.0, meaning no culling)

        std::vector<double> MolarAbundance() const;
    };

    struct NetOut {
        fourdst::composition::Composition composition; ///< Composition of the network after evaluation
        int num_steps; ///< Number of steps taken in the evaluation
        double energy; ///< Energy in ergs after evaluation

        friend std::ostream& operator<<(std::ostream& os, const NetOut& netOut) {
            os << "NetOut(composition=" << netOut.composition << ", num_steps=" << netOut.num_steps << ", energy=" << netOut.energy << ")";
            return os;
        }
    };

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
            virtual NetOut evaluate(const NetIn &netIn) = 0;

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



} // namespace nuclearNetwork
