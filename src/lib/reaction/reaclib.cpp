#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/atomic/species.h"

#include "gridfire/reaction/reaclib.h"
#include "gridfire/reaction/reactions_data.h"
#include "gridfire/types/types.h"
#include "gridfire/exceptions/error_reaction.h"

#include <stdexcept>
#include <sstream>
#include <vector>
#include <string>
#include <format>
#include <expected>

namespace {
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
        return {startIt, ritr.base()};
    }

    enum class ReactionParseError {
        MissingOpenParenthesis,
        MissingCloseParenthesis,
        ParenthesesOutOfOrder,
        MissingComma,
        BadFormat // Generic error (e.g., from an out_of_range exception)
    };

    std::string error_to_string(const ReactionParseError err) {
        switch (err) {
            case ReactionParseError::MissingOpenParenthesis:
                return "Missing '('";
            case ReactionParseError::MissingCloseParenthesis:
                return "Missing ')'";
            case ReactionParseError::ParenthesesOutOfOrder:
                return "')' found before '('";
            case ReactionParseError::MissingComma:
                return "Missing ',' within parentheses";
            case ReactionParseError::BadFormat:
                return "Bad format (substr error)";
        }
        return "Unknown error";
    }

    std::expected<std::string, ReactionParseError> reverse_pe_name(const std::string& forwardPEName) noexcept {
        try {
            size_t pos_open_paren = forwardPEName.find('(');
            size_t pos_close_paren = forwardPEName.find(')');

            if (pos_open_paren == std::string::npos) {
                return std::unexpected(ReactionParseError::MissingOpenParenthesis);
            }
            if (pos_close_paren == std::string::npos) {
                return std::unexpected(ReactionParseError::MissingCloseParenthesis);
            }
            if (pos_close_paren < pos_open_paren) {
                return std::unexpected(ReactionParseError::ParenthesesOutOfOrder);
            }

            std::string target = forwardPEName.substr(0, pos_open_paren);
            std::string result = forwardPEName.substr(pos_close_paren + 1);

            std::string inner_content = forwardPEName.substr(
                pos_open_paren + 1,
                pos_close_paren - pos_open_paren - 1
            );

            size_t pos_comma = inner_content.find(',');

            if (pos_comma == std::string::npos) {
                return std::unexpected(ReactionParseError::MissingComma);
            }

            std::string projectiles = inner_content.substr(0, pos_comma);

            std::string ejectiles = inner_content.substr(pos_comma + 1);

            std::ostringstream oss;
            oss << result << "(" << ejectiles << "," << projectiles << ")" << target;

            return oss.str();

        } catch (const std::out_of_range&) {
            return std::unexpected(ReactionParseError::BadFormat);
        }
    }
}
namespace gridfire::reaclib {
    static std::unique_ptr<reaction::ReactionSet> s_all_reaclib_reactions_ptr = nullptr;

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

        std::vector<std::unique_ptr<reaction::Reaction>> reaction_list;
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

            auto rpName_revNormalized = std::string(rpName_sv);
            if (reverse) {
                auto result = reverse_pe_name(std::string(rpName_sv));
                if (!result) {
                    std::string msg = std::format("Error reversing stored projectile-ejectile name for marked reverse reaction ({})", error_to_string(result.error()));
                    throw exceptions::ReactionParsingError(msg, rpName_revNormalized);
                }
                rpName_revNormalized = result.value();
            }

            // Construct the Reaction object. We use rpName for both the unique ID and the human-readable name.
            reaction_list.emplace_back(std::make_unique<reaction::ReaclibReaction>(
                rpName_revNormalized,
                rpName_revNormalized,
                chapter,
                reactants,
                products,
                qValue,
                label_sv,
                rate_coeffs,
                reverse
            ));
        }

        // The ReactionSet takes the vector of all individual reactions.
        const reaction::ReactionSet reaction_set(std::move(reaction_list));

        // The LogicalReactionSet groups reactions by their peName, which is what we want.
        s_all_reaclib_reactions_ptr = std::make_unique<reaction::ReactionSet>(reaction::packReactionSet(reaction_set));

        s_initialized = true;
    }


    // --- Public Interface Implementation ---

    const reaction::ReactionSet &get_all_reaclib_reactions() {
        // This ensures that the initialization happens only on the first call.
        if (!s_initialized) {
            initializeAllReaclibReactions();
        }
        if (s_all_reaclib_reactions_ptr == nullptr) {
            throw std::runtime_error("Reaclib reactions have not been initialized.");
        }
        return *s_all_reaclib_reactions_ptr;
    }

    bool reaction_is_weak(const reaction::Reaction& reaction) {
        const std::vector<fourdst::atomic::Species>& reactants = reaction.reactants();
        const std::vector<fourdst::atomic::Species>& products = reaction.products();

        if (reactants.size() != products.size()) {
            return false;
        }

        if (reactants.size() != 1 || products.size() != 1) {
            return false;
        }

        if (std::floor(reactants[0].a()) != std::floor(products[0].a())) {
            return false;
        }

        return true;
    }
} // namespace gridfire::reaclib