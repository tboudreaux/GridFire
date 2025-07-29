#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/species.h"

#include "gridfire/reaction/reaclib.h"
#include "gridfire/reaction/reactions_data.h"
#include "gridfire/network.h"

#include <stdexcept>
#include <sstream>
#include <vector>
#include <string>

std::string trim_whitespace(const std::string& str) {
    auto startIt = str.begin();
    const auto endIt   = str.end();

    while (startIt != endIt && std::isspace(static_cast<unsigned char>(*startIt))) {
        ++startIt;
    }
    if (startIt == endIt) {
        return "";
    }
    const auto ritr = std::find_if(str.rbegin(), std::string::const_reverse_iterator(startIt),
                             [](const unsigned char ch){ return !std::isspace(ch); });
    return std::string(startIt, ritr.base());
}

namespace gridfire::reaclib {
    static reaction::LogicalReactionSet* s_all_reaclib_reactions_ptr = nullptr;

    #pragma pack(push, 1)
    struct ReactionRecord {
        int32_t chapter;
        double qValue;
        double coeffs[7];
        bool reverse;
        char label[8];
        char rpName[64];
        char reactants_str[128];
        char products_str[128];
    };
    #pragma pack(pop)

    std::ostream& operator<<(std::ostream& os, const ReactionRecord& r) {
        os << "Chapter: " << r.chapter
           << ", Q-value: " << r.qValue
           << ", Coefficients: [" << r.coeffs[0] << ", " << r.coeffs[1] << ", "
           << r.coeffs[2] << ", " << r.coeffs[3] << ", " << r.coeffs[4] << ", "
           << r.coeffs[5] << ", " << r.coeffs[6] << "]"
           << ", Reverse: " << (r.reverse ? "true" : "false")
           << ", Label: '" << std::string(r.label, strnlen(r.label, sizeof(r.label))) << "'"
           << ", RP Name: '" << std::string(r.rpName, strnlen(r.rpName, sizeof(r.rpName))) << "'"
           << ", Reactants: '" << std::string(r.reactants_str, strnlen(r.reactants_str, sizeof(r.reactants_str))) << "'"
           << ", Products: '" << std::string(r.products_str, strnlen(r.products_str, sizeof(r.products_str))) << "'";
        return os;
    }

    static std::vector<fourdst::atomic::Species> parseSpeciesString(const std::string_view str) {
        std::vector<fourdst::atomic::Species> result;
        std::stringstream ss{std::string(str)};
        std::string name;

        while (ss >> name) {
            // Trim whitespace that might be left over from the fixed-width char arrays
            const auto trimmed_name = trim_whitespace(name);
            if (trimmed_name.empty()) continue;

            auto it = fourdst::atomic::species.find(trimmed_name);
            if (it != fourdst::atomic::species.end()) {
                result.push_back(it->second);
            } else {
                // If a species is not found, it's a critical data error.
                throw std::runtime_error("Unknown species in reaction data: " + std::string(trimmed_name));
            }
        }
        return result;
    }

    static void initializeAllReaclibReactions() {
        if (s_initialized) {
            return;
        }

        // Cast the raw byte data to our structured record format.
        const auto* records = reinterpret_cast<const ReactionRecord*>(raw_reactions_data);
        constexpr size_t num_reactions = raw_reactions_data_len / sizeof(ReactionRecord);

        std::vector<reaction::Reaction> reaction_list;
        reaction_list.reserve(num_reactions);

        for (size_t i = 0; i < num_reactions; ++i) {
            const auto&[chapter, qValue, coeffs, reverse, label, rpName, reactants_str, products_str] = records[i];

            // The char arrays from the binary are not guaranteed to be null-terminated
            // if the string fills the entire buffer. We create null-terminated string_views.
            const std::string_view label_sv(label, strnlen(label, sizeof(label)));
            const std::string_view rpName_sv(rpName, strnlen(rpName, sizeof(rpName)));
            const std::string_view reactants_sv(reactants_str, strnlen(reactants_str, sizeof(reactants_str)));
            const std::string_view products_sv(products_str, strnlen(products_str, sizeof(products_str)));

            auto reactants = parseSpeciesString(reactants_sv);
            auto products = parseSpeciesString(products_sv);

            const reaction::RateCoefficientSet rate_coeffs = {
                coeffs[0], coeffs[1], coeffs[2],
                coeffs[3], coeffs[4], coeffs[5],
                coeffs[6]
            };

            // Construct the Reaction object. We use rpName for both the unique ID and the human-readable name.
            reaction_list.emplace_back(
                rpName_sv,
                rpName_sv,
                chapter,
                reactants,
                products,
                qValue,
                label_sv,
                rate_coeffs,
                reverse
            );
        }

        // The ReactionSet takes the vector of all individual reactions.
        const reaction::ReactionSet reaction_set(std::move(reaction_list));

        // The LogicalReactionSet groups reactions by their peName, which is what we want.
        s_all_reaclib_reactions_ptr = new reaction::LogicalReactionSet(
            reaction::packReactionSetToLogicalReactionSet(reaction_set)
        );

        s_initialized = true;
    }


    // --- Public Interface Implementation ---

    const reaction::LogicalReactionSet& get_all_reactions() {
        // This ensures that the initialization happens only on the first call.
        if (!s_initialized) {
            initializeAllReaclibReactions();
        }
        if (s_all_reaclib_reactions_ptr == nullptr) {
            throw std::runtime_error("Reaclib reactions have not been initialized.");
        }
        return *s_all_reaclib_reactions_ptr;
    }
} // namespace gridfire::reaclib