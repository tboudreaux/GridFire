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

#include <ranges>

#include "gridfire/approx8.h"
#include "quill/LogMacros.h"
#include "gridfire/reaclib.h"
#include "gridfire/reactions.h"

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

    LogicalReaction::LogicalReaction(const std::vector<reaclib::REACLIBReaction> &reactions) {
        if (reactions.empty()) {
            LOG_ERROR(m_logger, "Empty reaction set provided to LogicalReaction constructor.");
            throw std::runtime_error("Empty reaction set provided to LogicalReaction constructor.");
        }

        const auto& first_reaction = reactions.front();
        m_peID = first_reaction.peName();
        m_reactants = first_reaction.reactants();
        m_products = first_reaction.products();
        m_qValue = first_reaction.qValue();
        m_reverse = first_reaction.is_reverse();
        m_chapter = first_reaction.chapter();
        m_rates.reserve(reactions.size());
        for (const auto& reaction : reactions) {
            m_rates.push_back(reaction.rateFits());
            if (std::abs(reaction.qValue() - m_qValue) > 1e-6) {
                LOG_ERROR(m_logger, "Inconsistent Q-values in reactions {}. All reactions must have the same Q-value.", m_peID);
                throw std::runtime_error("Inconsistent Q-values in reactions. All reactions must have the same Q-value.");
            }
        }
    }

    LogicalReaction::LogicalReaction(const reaclib::REACLIBReaction &first_reaction) {
        m_peID = first_reaction.peName();
        m_reactants = first_reaction.reactants();
        m_products = first_reaction.products();
        m_qValue = first_reaction.qValue();
        m_reverse = first_reaction.is_reverse();
        m_chapter = first_reaction.chapter();
        m_rates.reserve(1);
        m_rates.push_back(first_reaction.rateFits());

    }

    void LogicalReaction::add_reaction(const reaclib::REACLIBReaction &reaction) {
        if (std::abs(reaction.qValue() - m_qValue > 1e-6)) {
            LOG_ERROR(m_logger, "Inconsistent Q-values in reactions {}. All reactions must have the same Q-value.", m_peID);
            throw std::runtime_error("Inconsistent Q-values in reactions. All reactions must have the same Q-value.");
        }
        m_rates.push_back(reaction.rateFits());
    }

    auto LogicalReaction::begin() {
        return m_rates.begin();
    }
    auto LogicalReaction::begin() const {
        return m_rates.cbegin();
    }
    auto LogicalReaction::end() {
        return m_rates.end();
    }
    auto LogicalReaction::end() const {
        return m_rates.cend();
    }

    LogicalReactionSet::LogicalReactionSet(const std::vector<LogicalReaction> &reactions) {
        if (reactions.empty()) {
            LOG_ERROR(m_logger, "Empty reaction set provided to LogicalReactionSet constructor.");
            throw std::runtime_error("Empty reaction set provided to LogicalReactionSet constructor.");
        }

        for (const auto& reaction : reactions) {
            m_reactions.push_back(reaction);
        }

        m_reactionNameMap.reserve(m_reactions.size());
        for (const auto& reaction : m_reactions) {
            if (m_reactionNameMap.contains(reaction.id())) {
                LOG_ERROR(m_logger, "Duplicate reaction ID '{}' found in LogicalReactionSet.", reaction.id());
                throw std::runtime_error("Duplicate reaction ID found in LogicalReactionSet: " + std::string(reaction.id()));
            }
            m_reactionNameMap[reaction.id()] = reaction;
        }
    }

    LogicalReactionSet::LogicalReactionSet(const std::vector<reaclib::REACLIBReaction> &reactions) {
        std::vector<std::string_view> uniquePeNames;
        for (const auto& reaction: reactions) {
            if (uniquePeNames.empty()) {
                uniquePeNames.emplace_back(reaction.peName());
            } else if (std::ranges::find(uniquePeNames, reaction.peName()) == uniquePeNames.end()) {
                uniquePeNames.emplace_back(reaction.peName());
            }
        }

        for (const auto& peName : uniquePeNames) {
            std::vector<reaclib::REACLIBReaction> reactionsForPe;
            for (const auto& reaction : reactions) {
                if (reaction.peName() == peName) {
                    reactionsForPe.push_back(reaction);
                }
            }
            m_reactions.emplace_back(reactionsForPe);
        }

        m_reactionNameMap.reserve(m_reactions.size());
        for (const auto& reaction : m_reactions) {
            if (m_reactionNameMap.contains(reaction.id())) {
                LOG_ERROR(m_logger, "Duplicate reaction ID '{}' found in LogicalReactionSet.", reaction.id());
                throw std::runtime_error("Duplicate reaction ID found in LogicalReactionSet: " + std::string(reaction.id()));
            }
            m_reactionNameMap[reaction.id()] = reaction;
        }
    }

    LogicalReactionSet::LogicalReactionSet(const reaclib::REACLIBReactionSet &reactionSet) {
        LogicalReactionSet(reactionSet.get_reactions());
    }

    void LogicalReactionSet::add_reaction(const LogicalReaction& reaction) {
        if (m_reactionNameMap.contains(reaction.id())) {
            LOG_WARNING(m_logger, "Reaction {} already exists in the set. Not adding again.", reaction.id());
            std::cerr << "Warning: Reaction " << reaction.id() << " already exists in the set. Not adding again." << std::endl;
            return;
        }
        m_reactions.push_back(reaction);
        m_reactionNameMap[reaction.id()] = reaction;
    }

    void LogicalReactionSet::add_reaction(const reaclib::REACLIBReaction& reaction) {
        if (m_reactionNameMap.contains(reaction.id())) {
            m_reactionNameMap[reaction.id()].add_reaction(reaction);
        } else {
            const LogicalReaction logicalReaction(reaction);
            m_reactions.push_back(logicalReaction);
            m_reactionNameMap[reaction.id()] = logicalReaction;
        }
    }

    bool LogicalReactionSet::contains(const std::string_view &id) const {
        for (const auto& reaction : m_reactions) {
            if (reaction.id() == id) {
                return true;
            }
        }
        return false;
    }

    bool LogicalReactionSet::contains(const LogicalReaction &reactions) const {
        for (const auto& reaction : m_reactions) {
            if (reaction.id() == reactions.id()) {
                return true;
            }
        }
        return false;
    }

    bool LogicalReactionSet::contains_species(const fourdst::atomic::Species &species) const {
        return contains_reactant(species) || contains_product(species);
    }

    bool LogicalReactionSet::contains_reactant(const fourdst::atomic::Species &species) const {
        for (const auto& reaction : m_reactions) {
            if (std::ranges::find(reaction.reactants(), species) != reaction.reactants().end()) {
                return true;
            }
        }
        return false;
    }

    bool LogicalReactionSet::contains_product(const fourdst::atomic::Species &species) const {
        for (const auto& reaction : m_reactions) {
            if (std::ranges::find(reaction.products(), species) != reaction.products().end()) {
                return true;
            }
        }
        return false;
    }

    const LogicalReaction & LogicalReactionSet::operator[](size_t index) const {
        return m_reactions.at(index);
    }

    const LogicalReaction & LogicalReactionSet::operator[](const std::string_view &id) const {
        return m_reactionNameMap.at(id);
    }

    auto LogicalReactionSet::begin() {
        return m_reactions.begin();
    }

    auto LogicalReactionSet::begin() const {
        return m_reactions.cbegin();
    }

    auto LogicalReactionSet::end() {
        return m_reactions.end();
    }

    auto LogicalReactionSet::end() const {
        return m_reactions.cend();
    }

    LogicalReactionSet build_reaclib_nuclear_network(const fourdst::composition::Composition &composition) {
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
                LOG_INFO(logger, "Adding reaction {} to REACLIB reaction set.", reaction.peName());
                reactions.add_reaction(reaction);
            }
        }
        reactions.sort();
        return LogicalReactionSet(reactions);
    }

    LogicalReactionSet build_reaclib_nuclear_network(const fourdst::composition::Composition &composition, const double culling, const double T9) {
        using namespace reaclib;
        LogicalReactionSet allReactions = build_reaclib_nuclear_network(composition);
        LogicalReactionSet reactions;
        for (const auto& reaction : allReactions) {
            if (reaction.calculate_rate(T9) >= culling) {
                reactions.add_reaction(reaction);
            }
        }
        return reactions;
    }

}
