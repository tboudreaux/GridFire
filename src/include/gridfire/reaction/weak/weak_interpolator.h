#pragma once

#include "gridfire/reaction/weak/weak_types.h"
#include "fourdst/composition/atomicSpecies.h"

#include <unordered_map>
#include <cstdint>
#include <vector>
#include <expected>



namespace gridfire::rates::weak {
    class WeakRateInterpolator {
    public:
        using RowDataTable = std::array<RateDataRow, 77400>; // Total number of entries in the weak rate table NOTE: THIS MUST EQUAL THE VALUE IN weak_rate_library.h

        explicit WeakRateInterpolator(const RowDataTable& raw_data);

        [[nodiscard]] std::vector<fourdst::atomic::Species> available_isotopes() const;

        [[nodiscard]] std::expected<WeakRatePayload, InterpolationError> get_rates(
            uint16_t A,
            uint8_t Z,
            double t9,
            double log_rhoYe,
            double mu_e
        ) const;

        [[nodiscard]] std::expected<WeakRateDerivatives, InterpolationError> get_rate_derivatives(
            uint16_t A,
            uint8_t Z,
            double t9,
            double log_rhoYe,
            double mu_e
        ) const;
    private:
        static uint32_t pack_isotope_id(uint16_t A, uint8_t Z);

        std::unordered_map<uint32_t, IsotopeGrid> m_rate_table;
    };


}