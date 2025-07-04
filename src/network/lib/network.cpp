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
#include "gridfire/network.h"
#include "gridfire/reaction/reaclib.h"
#include "gridfire/reaction/reaction.h"

#include <ranges>

#include "quill/LogMacros.h"

namespace gridfire {
    std::vector<double> NetIn::MolarAbundance() const {
        std::vector <double> y;
        y.reserve(composition.getRegisteredSymbols().size());
        const auto [fst, snd] = composition.getComposition();
        for (const auto &name: fst | std::views::keys) {
            y.push_back(composition.getMolarAbundance(name));
        }
        return y;
    }


    Network::Network(const NetworkFormat format) :
        m_config(fourdst::config::Config::getInstance()),
        m_logManager(fourdst::logging::LogManager::getInstance()),
        m_logger(m_logManager.getLogger("log")),
        m_format(format),
        m_constants(fourdst::constant::Constants::getInstance()){
        if (format == NetworkFormat::UNKNOWN) {
            LOG_ERROR(m_logger, "nuclearNetwork::Network::Network() called with UNKNOWN format");
            m_logger->flush_log();
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

    reaction::LogicalReactionSet build_reaclib_nuclear_network(const fourdst::composition::Composition &composition, bool reverse) {
        using namespace reaction;
        std::vector<Reaction> reaclibReactions;
        auto logger = fourdst::logging::LogManager::getInstance().getLogger("log");

        for (const auto &reaction: reaclib::get_all_reactions()) {
            if (reaction.is_reverse() != reverse) {
                continue; // Skip reactions that do not match the requested direction
            }
            bool gotReaction = true;
            const auto& reactants = reaction.reactants();
            for (const auto& reactant : reactants) {
                if (!composition.contains(reactant)) {
                    gotReaction = false;
                    break; // If any reactant is not in the composition, skip this reaction
                }
            }
            if (gotReaction) {
                LOG_TRACE_L3(logger, "Adding reaction {} to REACLIB reaction set.", reaction.peName());
                reaclibReactions.push_back(reaction);
            }
        }
        const ReactionSet reactionSet(reaclibReactions);
        return packReactionSetToLogicalReactionSet(reactionSet);
    }

    // Trim whitespace from both ends of a string
    std::string trim_whitespace(const std::string& str) {
        auto startIt = str.begin();
        auto endIt   = str.end();

        while (startIt != endIt && std::isspace(static_cast<unsigned char>(*startIt))) {
            ++startIt;
        }
        if (startIt == endIt) {
            return "";
        }
        auto ritr = std::find_if(str.rbegin(), std::string::const_reverse_iterator(startIt),
                                 [](unsigned char ch){ return !std::isspace(ch); });
        return std::string(startIt, ritr.base());
    }

}
