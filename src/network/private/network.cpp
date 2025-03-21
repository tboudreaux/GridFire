#include "network.h"
#include "probe.h"
#include "quill/LogMacros.h"

namespace nuclearNetwork {
    Network::Network() :
        m_config(Config::getInstance()),
        m_logManager(Probe::LogManager::getInstance()),
        m_logger(m_logManager.getLogger("log")) {
    }
    nuclearNetwork::NetOut nuclearNetwork::Network::evaluate(const NetIn &netIn) {
        // You can throw an exception here or log a warning if it should never be used.
        LOG_ERROR(m_logger, "nuclearNetwork::Network::evaluate() is not implemented");
        throw std::runtime_error("nuclearNetwork::Network::evaluate() is not implemented");
    }
}
