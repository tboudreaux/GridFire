#pragma once

#include "gridfire/engine/engine_abstract.h"

#include <string>
#include <vector>

namespace gridfire::utils {
    /**
     * @brief Formats a map of nuclear species timescales into a human-readable string.
     *
     * This function takes a reaction network engine and the current plasma
     * conditions to calculate the characteristic timescales for each species.
     * It then formats this information into a neatly aligned ASCII table, which
     * is suitable for logging or printing to the console.
     *
     * @param engine A constant reference to a `DynamicEngine` object, used to
     *               calculate the species timescales.
     * @param Y A vector of the molar abundances (mol/g) for each species.
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
        const DynamicEngine& engine,
        const std::vector<double>& Y,
        const double T9,
        const double rho
    );
}