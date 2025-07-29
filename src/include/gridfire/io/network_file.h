#pragma once

#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"

#include "quill/Logger.h"

#include <string>
#include <vector>

namespace gridfire::io {
    typedef std::vector<std::string> ParsedNetworkData;

    /**
     * @class NetworkFileParser
     * @brief An abstract base class for network file parsers.
     *
     * This class defines the interface for parsing files that contain
     * reaction network definitions. Derived classes must implement the `parse`
     * method to handle specific file formats.
     */
    class NetworkFileParser {
    public:
        /**
         * @brief Virtual destructor for the base class.
         */
        virtual ~NetworkFileParser() = default;

        /**
         * @brief Parses a network file and returns the parsed data.
         *
         * This is a pure virtual function that must be implemented by derived
         * classes. It takes a filename as input and returns a `ParsedNetworkData`
         * struct containing the information extracted from the file.
         *
         * @param filename The path to the network file to parse.
         * @return A `ParsedNetworkData` struct containing the parsed reaction data.
         *
         * @throws std::runtime_error If the file cannot be opened or a parsing
         * error occurs.
         *
         * @b Usage
         * @code
         * std::unique_ptr<NetworkFileParser> parser = std::make_unique<SimpleReactionListFileParser>();
         * try {
         *     ParsedNetworkData data = parser->parse("my_reactions.txt");
         *     for (const auto& reaction_name : data.reactionPENames) {
         *         // ... process reaction name
         *     }
         * } catch (const std::runtime_error& e) {
         *     // ... handle error
         * }
         * @endcode
         */
        [[nodiscard]] virtual ParsedNetworkData parse(const std::string& filename) const = 0;

    };

    /**
     * @class SimpleReactionListFileParser
     * @brief A parser for simple text files containing a list of reactions.
     *
     * This parser reads a file where each line contains a single reaction name.
     * It supports comments (lines starting with '#') and ignores empty lines.
     *
     * @implements NetworkFileParser
     */
    class SimpleReactionListFileParser final : public NetworkFileParser {
    public:
        /**
         * @brief Constructs a SimpleReactionListFileParser.
         *
         * @post The parser is initialized and ready to parse files.
         */
        explicit SimpleReactionListFileParser();
        /**
         * @brief Parses a simple reaction list file.
         *
         * This method reads the specified file line by line. It trims whitespace
         * from each line, ignores lines that are empty or start with a '#'
         * comment character, and stores the remaining lines as reaction names.
         *
         * @param filename The path to the simple reaction list file.
         * @return A `ParsedNetworkData` struct containing the list of reaction names.
         *
         * @throws std::runtime_error If the file cannot be opened for reading.
         *
         * @b Algorithm
         * 1. Opens the specified file.
         * 2. Reads the file line by line.
         * 3. For each line, it removes any trailing comments (starting with '#').
         * 4. Trims leading and trailing whitespace.
         * 5. If the line is not empty, it is added to the list of reaction names.
         * 6. Returns the populated `ParsedNetworkData` struct.
         *
         * @b Usage
         * @code
         * SimpleReactionListFileParser parser;
         * ParsedNetworkData data = parser.parse("reactions.txt");
         * @endcode
         */
        ParsedNetworkData parse(const std::string& filename) const override;
    private:
        using Config = fourdst::config::Config;
        using LogManager = fourdst::logging::LogManager;
        Config& m_config = Config::getInstance();
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");
    };

    /**
     * @class MESANetworkFileParser
     * @brief A parser for MESA-format network files.
     *
     * This class is designed to parse reaction network files that follow the
     * format used by the MESA stellar evolution code.
     *
     * @implements NetworkFileParser
     */
    class MESANetworkFileParser final : public NetworkFileParser {
    public:
        /**
         * @brief Constructs a MESANetworkFileParser.
         *
         * @param filename The path to the MESA network file. This may be used
         * to pre-configure the parser.
         *
         * @post The parser is initialized with the context of the given file.
         */
        explicit MESANetworkFileParser(const std::string& filename);
        /**
         * @brief Parses a MESA-format network file.
         *
         * This method will read and interpret the structure of a MESA network
         * file to extract the list of reactions.
         *
         * @param filename The path to the MESA network file.
         * @return A `ParsedNetworkData` struct containing the list of reaction names.
         *
         * @throws std::runtime_error If the file cannot be opened or if it
         * contains formatting errors.
         */
        ParsedNetworkData parse(const std::string& filename) const override;
    private:
        using Config = fourdst::config::Config;
        using LogManager = fourdst::logging::LogManager;
        Config& m_config = Config::getInstance();
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");

        std::string m_filename;
    };
}