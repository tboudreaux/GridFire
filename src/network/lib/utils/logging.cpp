#include "gridfire/utils/logging.h"
#include "gridfire/engine/engine_abstract.h"

#include <sstream>
#include <iomanip>
#include <algorithm>
#include <ranges>
#include <string_view>
#include <string>
#include <iostream>
#include <vector>

std::string gridfire::utils::formatNuclearTimescaleLogString(
    const DynamicEngine& engine,
    std::vector<double> const& Y,
    const double T9,
    const double rho
) {
    auto const& result = engine.getSpeciesTimescales(Y, T9, rho);
    if (!result) {
        std::ostringstream ss;
        ss << "Failed to get species timescales: " << result.error();
        return ss.str();
    }
    const std::unordered_map<fourdst::atomic::Species, double>& timescales = result.value();

    // Figure out how wide the "Species" column needs to be:
    std::size_t maxNameLen = std::string_view("Species").size();
    for (const auto &key: timescales | std::views::keys) {
        std::string_view name = key.name();
        maxNameLen = std::max(maxNameLen, name.size());
    }

    // Pick a fixed width for the timescale column:
    constexpr int timescaleWidth = 12;

    std::ostringstream ss;
    ss << "== Timescales (s) ==\n";

    // Header row
    ss << std::left << std::setw(static_cast<int>(maxNameLen) + 2) << "Species"
       << std::right << std::setw(timescaleWidth)     << "Timescale (s)" << "\n";

    // Underline
    ss << std::string(static_cast<int>(maxNameLen) + 2 + timescaleWidth, '=') << "\n";

    ss << std::scientific;

    // Data rows
    for (auto const& [species, timescale] : timescales) {
        const std::string_view name = species.name();
        ss << std::left << std::setw(static_cast<int>(maxNameLen) + 2) << name;

        if (std::isinf(timescale)) {
            ss << std::right << std::setw(timescaleWidth) << "inf" << "\n";
        } else {
            ss << std::right << std::setw(timescaleWidth)
               << std::scientific << std::setprecision(3) << timescale << "\n";
        }
    }

    // Footer underline
    ss << std::string(static_cast<int>(maxNameLen) + 2 + timescaleWidth, '=') << "\n";

    return ss.str();
}
