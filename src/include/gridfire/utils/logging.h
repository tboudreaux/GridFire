#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "fourdst/composition/composition.h"
#include "gridfire/engine/scratchpads/blob.h"

#include <string>
#include <functional>

namespace gridfire::utils {
    /**
     * @brief Formats a map of nuclear species timescales into a human-readable string.
     *
     * This function takes a reaction network engine and the current plasma
     * conditions to calculate the characteristic timescales for each species.
     * It then formats this information into a neatly aligned ASCII table, which
     * is suitable for logging or printing to the console.
     *
     * @param ctx
     * @param engine A constant reference to a `DynamicEngine` object, used to
     *               calculate the species timescales.
     * @param composition The current composition of the plasma
     * @param T9 The temperature in units of 10^9 K.
     * @param rho The plasma density in g/cm^3.
     * @return A std::string containing the formatted table of species and their
     *         timescales.
     *
     * @b Pre-conditions
     * - The `engine` must be in a valid state.
     * - The size of the `Y` vector must be consistent with the number of species
     *   expected by the `engine`.
     *
     * @b Algorithm
     * 1. Calls the `getSpeciesTimescales` method on the provided `engine` to get
     *    the timescale for each species under the given conditions.
     * 2. Determines the maximum length of the species names to dynamically set the
     *    width of the "Species" column for proper alignment.
     * 3. Uses a `std::ostringstream` to build the output string.
     * 4. Constructs a header for the table with titles "Species" and "Timescale (s)".
     * 5. Iterates through the map of timescales, adding a row to the table for
     *    each species.
     * 6. Timescales are formatted in scientific notation with 3 digits of precision.
     * 7. Special handling is included to print "inf" for infinite timescales.
     * 8. The final string, including header and footer lines, is returned.
     *
     * @b Usage
     * @code
     * // Assume 'my_engine' is a valid DynamicEngine object and Y, T9, rho are initialized.
     * std::string log_output = gridfire::utils::formatNuclearTimescaleLogString(my_engine, Y, T9, rho);
     * std::cout << log_output;
     *
     * // Example Output:
     * // == Timescales (s) ==
     * // Species      Timescale (s)
     * // ==========================
     * // h1           1.234e+05
     * // he4          inf
     * // c12          8.765e-02
     * // ==========================
     * @endcode
     */
    std::string formatNuclearTimescaleLogString(
        engine::scratch::StateBlob &ctx,
        const engine::DynamicEngine& engine,
        const fourdst::composition::Composition& composition,
        double T9, double rho
    );

    template <typename T>
    concept Streamable = requires(std::ostream& os, const T& value) {
            { os << value } -> std::same_as<std::ostream&>;
    };

    template <
        std::ranges::input_range Container,
        typename Elem = std::ranges::range_reference_t<Container>,
        typename Transform = std::identity,
        typename Pred = bool(*)(const std::ranges::range_value_t<Container>&)
    >
    requires std::invocable<Transform, Elem> && Streamable<std::invoke_result_t<Transform, Elem>> && std::predicate<Pred, Elem>
    static std::string iterable_to_delimited_string(
        const Container& container,
        const std::string_view delimiter = ", ",
        Transform transform = {},
        Pred pred = [](const auto&){ return true; }
    ) noexcept {
        std::ostringstream oss;
        bool first = true;
        for (auto&& item : container) {
            if (!std::invoke(pred, item)) {
                continue;
            }
            if (!first) {
                oss << delimiter;
            }
            oss << std::invoke(transform, item);
            first = false;
        }
        return oss.str();
    }




}