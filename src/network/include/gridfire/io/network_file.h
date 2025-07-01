#pragma once

#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"

#include "quill/Logger.h"

#include <string>
#include <vector>

namespace gridfire::io {

    struct ParsedNetworkData {
        std::vector<std::string> reactionPENames;
    };

    class NetworkFileParser {
    public:
        virtual ~NetworkFileParser() = default;

        [[nodiscard]] virtual ParsedNetworkData parse(const std::string& filename) const = 0;

    };

    class SimpleReactionListFileParser final : public NetworkFileParser {
    public:
        explicit SimpleReactionListFileParser();
        ParsedNetworkData parse(const std::string& filename) const override;
    private:
        using Config = fourdst::config::Config;
        using LogManager = fourdst::logging::LogManager;
        Config& m_config = Config::getInstance();
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");
    };

    class MESANetworkFileParser final : public NetworkFileParser {
    public:
        explicit MESANetworkFileParser(const std::string& filename);
        ParsedNetworkData parse(const std::string& filename) const override;
    private:
        using Config = fourdst::config::Config;
        using LogManager = fourdst::logging::LogManager;
        Config& m_config = Config::getInstance();
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");

        std::string m_filename;
    };
}