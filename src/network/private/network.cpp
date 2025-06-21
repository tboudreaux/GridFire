/* ***********************************************************************
//
//   Copyright (C) 2025 -- The 4D-STAR Collaboration
//   File Authors: Aaron Dotter, Emily Boudreaux
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
#include "network.h"

#include <ranges>

#include "approx8.h"
#include "quill/LogMacros.h"
#include "reaclib.h"
#include "reactions.h"

namespace gridfire {
    Network::Network(const NetworkFormat format) :
        m_config(fourdst::config::Config::getInstance()),
        m_logManager(fourdst::logging::LogManager::getInstance()),
        m_logger(m_logManager.getLogger("log")),
        m_format(format),
        m_constants(fourdst::constant::Constants::getInstance()){
        if (format == NetworkFormat::UNKNOWN) {
            LOG_ERROR(m_logger, "nuclearNetwork::Network::Network() called with UNKNOWN format");
            throw std::runtime_error("nuclearNetwork::Network::Network() called with UNKNOWN format");
        }
    }

    NetworkFormat Network::getFormat() const {
        return m_format;
    }

    NetworkFormat Network::setFormat(const NetworkFormat format) {
        const NetworkFormat oldFormat = m_format;
        m_format = format;
        return oldFormat;
    }

    NetOut Network::evaluate(const NetIn &netIn) {
        NetOut netOut;
        switch (m_format) {
            case APPROX8: {
                approx8::Approx8Network network;
                netOut = network.evaluate(netIn);
                break;
            }
            case UNKNOWN: {
                LOG_ERROR(m_logger, "Network format {} is not implemented.", FormatStringLookup.at(m_format));
                throw std::runtime_error("Network format not implemented.");
            }
            default: {
                LOG_ERROR(m_logger, "Unknown network format.");
                throw std::runtime_error("Unknown network format.");
            }
        }
        return netOut;
    }

    reaclib::REACLIBReactionSet build_reaclib_nuclear_network(const fourdst::composition::Composition &composition) {
        using namespace reaclib;
        REACLIBReactionSet reactions;
        auto logger = fourdst::logging::LogManager::getInstance().getLogger("log");

        if (!s_initialized) {
            LOG_INFO(logger, "REACLIB reactions not initialized. Calling initializeAllReaclibReactions()...");
            initializeAllReaclibReactions();
        }

        for (const auto &reaction: s_all_reaclib_reactions | std::views::values) {
            bool gotReaction = true;
            const auto& reactants = reaction.reactants();
            for (const auto& reactant : reactants) {
                if (!composition.contains(reactant)) {
                    gotReaction = false;
                    break; // If any reactant is not in the composition, skip this reaction
                }
            }
            if (gotReaction) {
                LOG_INFO(logger, "Adding reaction {} to REACLIB reaction set.", reaction.id());
                reactions.add_reaction(reaction);
            }
        }
        reactions.sort();
        return reactions;
    }

    reaclib::REACLIBReactionSet build_reaclib_nuclear_network(const fourdst::composition::Composition &composition, const double culling, const double T9) {
        using namespace reaclib;
        REACLIBReactionSet allReactions = build_reaclib_nuclear_network(composition);
        REACLIBReactionSet reactions;
        for (const auto& reaction : allReactions) {
            if (reaction.calculate_rate(T9) >= culling) {
                reactions.add_reaction(reaction);
            }
        }
        return reactions;
    }

}
