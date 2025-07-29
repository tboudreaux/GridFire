#include "gridfire/io/network_file.h"

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <stdexcept>

#include "quill/LogMacros.h"

namespace gridfire::io {
    namespace {
        inline void ltrim(std::string &s) {
           s.erase(
               s.begin(),
               std::ranges::find_if(s,
                                    [](const unsigned char ch) {
                                        return !std::isspace(ch);
                                    })
            );
        }

        inline void rtrim(std::string &s) {
            s.erase(
                std::find_if(
                    s.rbegin(),
                    s.rend(),
                    [](const unsigned char ch) {
                        return !std::isspace(ch);
                    }).base(),
                    s.end()
                );
        }

        inline void trim(std::string &s) {
            ltrim(s);
            rtrim(s);
        }


    }
    SimpleReactionListFileParser::SimpleReactionListFileParser() {}

    ParsedNetworkData SimpleReactionListFileParser::parse(const std::string& filename) const {
        LOG_TRACE_L1(m_logger, "Parsing simple reaction list file: {}", filename);

        std::ifstream file(filename);
        if (!file.is_open()) {
            LOG_ERROR(m_logger, "Failed to open file: {}", filename);
            m_logger -> flush_log();
            throw std::runtime_error("Could not open file: " + filename);
        }

        ParsedNetworkData parsed;
        std::string line;
        while (std::getline(file, line)) {
            LOG_TRACE_L3(m_logger, "Parsing reaction list file {}, line: {}", filename, line);

            const size_t comment_pos = line.find('#');
            if (comment_pos != std::string::npos) {
                line = line.substr(0, comment_pos);
            }

            trim(line);

            if (line.empty()) {
                continue; // Skip empty lines
            }
            parsed.push_back(line);
        }
        LOG_TRACE_L1(m_logger, "Parsed {} reactions from file: {}", parsed.size(), filename);
        return parsed;
    }

}