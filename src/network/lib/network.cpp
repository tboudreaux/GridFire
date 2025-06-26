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
#include "gridfire/reactions.h"
#include "../include/gridfire/reaction/reaction.h"

#include <ranges>
#include <fstream>

#include "quill/LogMacros.h"

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

    reaction::REACLIBLogicalReactionSet build_reaclib_nuclear_network(const fourdst::composition::Composition &composition, bool reverse) {
        using namespace reaction;
        std::vector<reaction::REACLIBReaction> reaclibReactions;
        auto logger = fourdst::logging::LogManager::getInstance().getLogger("log");

        if (!reaclib::s_initialized) {
            LOG_DEBUG(logger, "REACLIB reactions not initialized. Calling initializeAllReaclibReactions()...");
            reaclib::initializeAllReaclibReactions();
        }

        for (const auto &reaction: reaclib::s_all_reaclib_reactions | std::views::values) {
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
        const REACLIBReactionSet reactionSet(reaclibReactions);
        return REACLIBLogicalReactionSet(reactionSet);
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

    reaction::REACLIBLogicalReactionSet build_reaclib_nuclear_network_from_file(const std::string& filename, bool reverse) {
        const auto logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::vector<std::string> reactionPENames;
        std::ifstream infile(filename, std::ios::in);
        if (!infile.is_open()) {
            LOG_ERROR(logger, "Failed to open network definition file: {}", filename);
            throw std::runtime_error("Failed to open network definition file: " + filename);
        }
        std::string line;
        while (std::getline(infile, line)) {
            std::string trimmedLine = trim_whitespace(line);
            if (trimmedLine.empty() || trimmedLine.starts_with('#')) {
                continue; // Skip empty lines and comments
            }
            reactionPENames.push_back(trimmedLine);
        }
        infile.close();
        std::vector<reaction::REACLIBReaction> reaclibReactions;

        if (reactionPENames.empty()) {
            LOG_ERROR(logger, "No reactions found in the network definition file: {}", filename);
            throw std::runtime_error("No reactions found in the network definition file: " + filename);
        }

        if (!reaclib::s_initialized) {
            LOG_DEBUG(logger, "REACLIB reactions not initialized. Calling initializeAllReaclibReactions()...");
            reaclib::initializeAllReaclibReactions();
        } else {
            LOG_DEBUG(logger, "REACLIB reactions already initialized.");
        }

        for (const auto& peName : reactionPENames) {
            bool found = false;
            for (const auto& reaction : reaclib::s_all_reaclib_reactions | std::views::values) {
                if (reaction.peName() == peName && reaction.is_reverse() == reverse) {
                    reaclibReactions.push_back(reaction);
                    found = true;
                    LOG_TRACE_L3(logger, "Found reaction {} (version {}) in REACLIB database.", peName, reaction.sourceLabel());
                }
            }
            if (!found) {
                LOG_ERROR(logger, "Reaction {} not found in REACLIB database. Skipping...", peName);
                throw std::runtime_error("Reaction not found in REACLIB database.");
            }
        }


        // Log any reactions that were not found in the REACLIB database
        for (const auto& peName : reactionPENames) {
            if (std::ranges::find(reaclibReactions, peName, &reaction::REACLIBReaction::peName) == reaclibReactions.end()) {
                LOG_WARNING(logger, "Reaction {} not found in REACLIB database.", peName);
            }
        }
        const reaction::REACLIBReactionSet reactionSet(reaclibReactions);
        return reaction::REACLIBLogicalReactionSet(reactionSet);
    }

}
