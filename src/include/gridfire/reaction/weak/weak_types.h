#pragma once

/**
 * @file weak_types.h
 * @brief Plain data structures and enums for weak reaction tables, interpolation payloads, and errors.
 *
 * This header defines the raw row format loaded from the unified weak-rate library, simple
 * enumerations for channels and axes, compact payloads for interpolated values and derivatives,
 * and error-reporting structures used by the interpolator.
 */

#include <cstdint>
#include <array>
#include <vector>
#include <optional>
#include <unordered_map>
#include <ostream>

namespace gridfire::rates::weak {
    /**
     * @brief One row of the unified weak-rate data table for a specific isotope and state.
     *
     * Units and meanings:
     *  - t9: temperature in GK (10^9 K).
     *  - log_rhoye: base-10 logarithm of rho*Ye where rho is g cm^-3 and Ye is electron fraction.
     *  - mu_e: electron chemical potential in MeV.
     *  - log_*: base-10 logarithm of the tabulated rate or neutrino-energy loss term.
     *
     * Channel mappings:
     *  - beta-plus (β+): log_beta_plus, neutrino-loss column log_neutrino_loss_ec.
     *  - electron capture (e− cap): log_electron_capture, neutrino-loss column log_neutrino_loss_ec.
     *  - beta-minus (β−): log_beta_minus, neutrino-loss column log_antineutrino_loss_bd.
     *  - positron capture (e+ cap): log_positron_capture, neutrino-loss column log_antineutrino_loss_bd.
     */
    struct RateDataRow {
        uint16_t A;   ///< Mass number.
        uint8_t Z;    ///< Proton number.
        float t9;     ///< Temperature in GK.
        float log_rhoye; ///< log10(rho*Ye) (cgs density times electron fraction).
        float mu_e;   ///< Electron chemical potential (MeV).
        float log_beta_plus;            ///< log10(β+ decay rate).
        float log_electron_capture;     ///< log10(e− capture rate).
        float log_neutrino_loss_ec;     ///< log10(neutrino loss for β+ and e− capture).
        float log_beta_minus;           ///< log10(β− decay rate).
        float log_positron_capture;     ///< log10(e+ capture rate).
        float log_antineutrino_loss_bd; ///< log10(antineutrino loss for β− and e+ capture).
    };

    /**
     * @brief Weak reaction channel identifiers.
     */
    enum class WeakReactionType {
        BETA_PLUS_DECAY,  ///< β+ decay: Z -> Z-1 + e+ + ν_e
        BETA_MINUS_DECAY, ///< β− decay: Z -> Z+1 + e− + ν̄_e
        ELECTRON_CAPTURE, ///< e− capture: (Z, e−) -> Z-1 + ν_e
        POSITRON_CAPTURE, ///< e+ capture: (Z, e+) -> Z+1 + ν̄_e
    };

    /**
     * @brief Enumeration of neutrino flavors (for potential extensions and tagging).
     */
    enum class NeutrinoTypes {
        ELECTRON_NEUTRINO,
        ELECTRON_ANTINEUTRINO,
        MUON_NEUTRINO,
        MUON_ANTINEUTRINO,
        TAU_NEUTRINO,
        TAU_ANTINEUTRINO
    };

    /**
     * @brief Lookup errors for WeakReactionMap queries.
     */
    enum class WeakMapError {
        SPECIES_NOT_FOUND, ///< No entries for the requested Species.
        UNKNOWN_ERROR
    };

    /**
     * @brief Interpolated weak-rate payload at a single state.
     *
     * All values are base-10 logarithms of the corresponding rates or neutrino-loss terms.
     * Consumers typically convert with pow(10, log_value) and may apply sentinel thresholds
     * at the usage site.
     */
    struct WeakRatePayload {
        double log_beta_plus;            ///< log10(β+ decay rate).
        double log_electron_capture;     ///< log10(e− capture rate).
        double log_neutrino_loss_ec;     ///< log10(neutrino loss for β+ and e− capture).
        double log_beta_minus;           ///< log10(β− decay rate).
        double log_positron_capture;     ///< log10(e+ capture rate).
        double log_antineutrino_loss_bd; ///< log10(antineutrino loss for β− and e+ capture).
    };

    /**
     * @brief Partial derivatives of the log10() fields w.r.t. (T9, log10(rho*Ye)).
     *
     * Array ordering is [d/dT9, d/dlogRhoYe] for each corresponding field.
     */
    struct WeakRateDerivatives {
        // Each array holds [d/dT9, d/dlogRhoYe, d/dMuE]
        std::array<double, 2> d_log_beta_plus;
        std::array<double, 2> d_log_electron_capture;
        std::array<double, 2> d_log_neutrino_loss_ec;
        std::array<double, 2> d_log_beta_minus;
        std::array<double, 2> d_log_positron_capture;
        std::array<double, 2> d_log_antineutrino_loss_bd;
    };

    /**
     * @brief Error categories for interpolation attempts.
     */
    enum class InterpolationErrorType {
        BOUNDS_ERROR,          ///< Query outside the per-axis min/max of the table.
        UNKNOWN_SPECIES_ERROR, ///< Requested (A,Z) not present in the tables.
        UNKNOWN_ERROR
    };

    /**
     * @brief Human-readable names for InterpolationErrorType.
     */
    inline std::unordered_map<InterpolationErrorType, std::string_view> InterpolationErrorTypeMap = {
        {InterpolationErrorType::BOUNDS_ERROR, "Bounds Error"},
        {InterpolationErrorType::UNKNOWN_SPECIES_ERROR, "Unknown Species Error"},
        {InterpolationErrorType::UNKNOWN_ERROR, "Unknown Error"}
    };

    /**
     * @brief Axes of the interpolation table.
     */
    enum class TableAxes {
        T9,        ///< Temperature in GK.
        LOG_RHOYE, ///< log10(rho*Ye).
        MUE        ///< Electron chemical potential (MeV).
    };
}

// This need to be here to avoid compiler issues related to the order of specialization
namespace std {
    template <>
    struct hash<gridfire::rates::weak::TableAxes> {
        std::size_t operator()(gridfire::rates::weak::TableAxes t) const noexcept {
            return std::hash<int>()(static_cast<int>(t));
        }
    };
}

namespace gridfire::rates::weak {

    /**
     * @brief Detailed bounds information for a BOUNDS_ERROR.
     */
    struct BoundsErrorInfo {
        TableAxes axis;      ///< Axis on which the error occurred.
        double axisMinValue; ///< Minimum tabulated value on the axis.
        double axisMaxValue; ///< Maximum tabulated value on the axis.
        double queryValue;   ///< Requested value.
    };

    /**
     * @brief Interpolation error with optional per-axis bounds details.
     *
     * For BOUNDS_ERROR, boundsErrorInfo may contain an entry per offending axis.
     */
    struct InterpolationError {
        InterpolationErrorType type; ///< Error category.
        std::optional<std::unordered_map<TableAxes, BoundsErrorInfo>> boundsErrorInfo = std::nullopt;
    };


    /**
     * @brief Regular 2D grid and payloads for a single isotope (A,Z).
     *
     * Axes are monotonically increasing per dimension. Data vector is laid out in
     * row-major order with index computed as:
     *
     *   index = i_t9 * N_rhoYe + j_rhoYe
     *
     */
    struct IsotopeGrid {
        std::vector<double> t9_axis;    ///< Unique sorted T9 grid.
        std::vector<double> rhoYe_axis; ///< Unique sorted log10(rho*Ye) grid.
        std::vector<WeakRatePayload> data; ///< MuE axis for each (T9, log_rhoYe) pair (the table is ragged in mu_e). This is also where the payloads are stored.
    };

    /**
     * @brief Abbreviated channel name used in printing and IDs.
     * @param t Channel enum.
     * @return Short name: bp, bm, ec, or pc.
     */
    constexpr std::string_view weak_reaction_type_name(const WeakReactionType t) noexcept {
        switch (t) {
            case WeakReactionType::BETA_PLUS_DECAY:   return "bp";
            case WeakReactionType::BETA_MINUS_DECAY:  return "bm";
            case WeakReactionType::ELECTRON_CAPTURE:  return "ec";
            case WeakReactionType::POSITRON_CAPTURE:  return "pc";
        }
        return "Unknown";
    }

    /**
     * @brief A single weak-reaction data point (type, state, and log values).
     *
     * All rates and losses are base-10 logarithms. Useful for listing and filtering
     * weak entries for a Species.
     *
     * @par Example
     * @code
     * WeakReactionEntry e{WeakReactionType::ELECTRON_CAPTURE, 3.0f, 6.0f, 2.0f, -2.3f, -1.7f};
     * std::cout << e << "\n"; // prints a compact summary
     * @endcode
     */
    struct WeakReactionEntry {
        WeakReactionType type; ///< Channel.
        float T9;              ///< Temperature in GK.
        float log_rhoYe;       ///< log10(rho*Ye).
        float mu_e;            ///< Electron chemical potential (MeV).
        float log_rate;        ///< Channel-specific log10(rate).
        float log_neutrino_loss; ///< Corresponding log10(neutrino or antineutrino energy loss).

        friend std::ostream& operator<<(std::ostream& os, const WeakReactionEntry& reaction) {
            os << "WeakReactionEntry(type=" << weak_reaction_type_name(reaction.type)
               << ", T9=" << reaction.T9
               << ", log_rhoYe=" << reaction.log_rhoYe
               << ", mu_e=" << reaction.mu_e
               << ", log_rate=" << reaction.log_rate
               << ", log_neutrino_loss=" << reaction.log_neutrino_loss
               << ")";
            return os;
        }
    };
}

