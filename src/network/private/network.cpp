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
#include "probe.h"
#include "quill/LogMacros.h"

namespace serif::network {
    Network::Network() :
        m_config(serif::config::Config::getInstance()),
        m_logManager(serif::probe::LogManager::getInstance()),
        m_logger(m_logManager.getLogger("log")) {
    }
    NetOut Network::evaluate(const NetIn &netIn) {
        // You can throw an exception here or log a warning if it should never be used.
        LOG_ERROR(m_logger, "nuclearNetwork::Network::evaluate() is not implemented");
        throw std::runtime_error("nuclearNetwork::Network::evaluate() is not implemented");
    }
}
