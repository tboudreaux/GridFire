#pragma once

#include <cstdint>
#include <array>
#include <vector>
#include <optional>
#include <unordered_map>
#include <ostream>

namespace gridfire::rates::weak {
    struct RateDataRow {
        uint16_t A;
        uint8_t Z;
        float t9;
        float log_rhoye;
        float mu_e;
        float log_beta_plus;
        float log_electron_capture;
        float log_neutrino_loss_ec;
        float log_beta_minus;
        float log_positron_capture;
        float log_antineutrino_loss_bd;
    };

    enum class WeakReactionType {
        BETA_PLUS_DECAY,
        BETA_MINUS_DECAY,
        ELECTRON_CAPTURE,
        POSITRON_CAPTURE,
    };

    enum class NeutrinoTypes {
        ELECTRON_NEUTRINO,
        ELECTRON_ANTINEUTRINO,
        MUON_NEUTRINO,
        MUON_ANTINEUTRINO,
        TAU_NEUTRINO,
        TAU_ANTINEUTRINO
    };

    enum class WeakMapError {
        SPECIES_NOT_FOUND,
        UNKNOWN_ERROR
    };

    struct WeakRatePayload {
        double log_beta_plus;
        double log_electron_capture;
        double log_neutrino_loss_ec;
        double log_beta_minus;
        double log_positron_capture;
        double log_antineutrino_loss_bd;
    };

    struct WeakRateDerivatives {
        // Each array holds [d/dT9, d/dlogRhoYe, d/dMuE]
        std::array<double, 3> d_log_beta_plus;
        std::array<double, 3> d_log_electron_capture;
        std::array<double, 3> d_log_neutrino_loss_ec;
        std::array<double, 3> d_log_beta_minus;
        std::array<double, 3> d_log_positron_capture;
        std::array<double, 3> d_log_antineutrino_loss_bd;
    };

    enum class InterpolationErrorType {
        BOUNDS_ERROR,
        UNKNOWN_SPECIES_ERROR,
        UNKNOWN_ERROR
    };

    inline std::unordered_map<InterpolationErrorType, std::string_view> InterpolationErrorTypeMap = {
        {InterpolationErrorType::BOUNDS_ERROR, "Bounds Error"},
        {InterpolationErrorType::UNKNOWN_SPECIES_ERROR, "Unknown Species Error"},
        {InterpolationErrorType::UNKNOWN_ERROR, "Unknown Error"}
    };

    enum class TableAxes {
        T9,
        LOG_RHOYE,
        MUE
    };

    struct BoundsErrorInfo {
        TableAxes axis;
        double axisMinValue;
        double axisMaxValue;
        double queryValue;
    };

    struct InterpolationError {
        InterpolationErrorType type;
        std::optional<std::unordered_map<TableAxes, BoundsErrorInfo>> boundsErrorInfo = std::nullopt;
    };

    struct IsotopeGrid {
        std::vector<double> t9_axis;
        std::vector<double> rhoYe_axis;
        std::vector<double> mue_axis;

        // index = (i_t9 * logRhoYe_axis.size() + j_rhoYe) + mue_axis.size() + k_mue
        std::vector<WeakRatePayload> data;
    };

    constexpr std::string_view weak_reaction_type_name(const WeakReactionType t) noexcept {
        switch (t) {
            case WeakReactionType::BETA_PLUS_DECAY:   return "bp";
            case WeakReactionType::BETA_MINUS_DECAY:  return "bm";
            case WeakReactionType::ELECTRON_CAPTURE:  return "ec";
            case WeakReactionType::POSITRON_CAPTURE:  return "pc";
        }
        return "Unknown";
    }

    struct WeakReactionEntry {
        WeakReactionType type;
        float T9;
        float log_rhoYe;
        float mu_e;
        float log_rate;
        float log_neutrino_loss;

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