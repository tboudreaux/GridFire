#include "gridfire/partition/partition_rauscher_thielemann.h"
#include "gridfire/partition/rauscher_thielemann_partition_data.h"
#include "gridfire/partition/rauscher_thielemann_partition_data_record.h"

#include "fourdst/logging/logging.h"
#include "quill/LogMacros.h"

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <array>
#include <iostream>

namespace gridfire::partition {
    static constexpr std::array<double, 24> RT_TEMPERATURE_GRID_T9 = {
        0.01, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5,
        2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0
    };

    RauscherThielemannPartitionFunction::RauscherThielemannPartitionFunction() {
        constexpr size_t numRecords = rauscher_thielemann_partition_data_len / sizeof(record::RauscherThielemannPartitionDataRecord);
        m_partitionData.reserve(numRecords);
        const auto* records = reinterpret_cast<const record::RauscherThielemannPartitionDataRecord*>(rauscher_thielemann_partition_data);
        for (size_t i = 0; i < numRecords; ++i) {
            const auto&[z, a, ground_state_spin, normalized_g_values] = records[i];
            IsotopeData data;
            data.ground_state_spin = ground_state_spin;
            std::ranges::copy(normalized_g_values, data.normalized_g_values.begin());
            const int key = make_key(z, a);
            LOG_TRACE_L3_LIMIT_EVERY_N(
                100,
                m_logger,
                "(EVERY 100) Adding Rauscher-Thielemann partition data for Z={} A={} (key={})",
                z,
                a,
                key
            );

            m_partitionData[key] = std::move(data);
        }
    }

    double RauscherThielemannPartitionFunction::evaluate(
        const int z,
        const int a,
        const double T9
    ) const {
        LOG_TRACE_L2(m_logger, "Evaluating Rauscher-Thielemann partition function for Z={} A={} T9={}", z, a, T9);

        const auto [bound, data, upperIndex, lowerIndex] = find(z, a, T9);

        switch (bound) {
            case FRONT: {
                LOG_TRACE_L2(m_logger, "Using FRONT bound for Z={} A={} T9={}", z, a, T9);
                return data.normalized_g_values.front() * (2.0 * data.ground_state_spin + 1.0);
            }
            case BACK: {
                LOG_TRACE_L2(m_logger, "Using BACK bound for Z={} A={} T9={}", z, a, T9);
                return data.normalized_g_values.back() * (2.0 * data.ground_state_spin + 1.0);
            }
            case MIDDLE: {
                LOG_TRACE_L2(m_logger, "Using MIDDLE bound for Z={} A={} T9={}", z, a, T9);
            }
        }

        const auto [T9_high, G_norm_high, T9_low, G_norm_low] = get_interpolation_points(
            upperIndex,
            lowerIndex,
            data.normalized_g_values
        );

        const double frac = (T9 - T9_low) / (T9_high - T9_low);
        const double interpolated_g_norm = G_norm_low + frac * (G_norm_high - G_norm_low);

        return interpolated_g_norm * (2.0 * data.ground_state_spin + 1.0);
    }

    double RauscherThielemannPartitionFunction::evaluateDerivative(
        const int z,
        const int a,
        const double T9
    ) const {
        LOG_TRACE_L2(m_logger, "Evaluating derivative of Rauscher-Thielemann partition function for Z={} A={} T9={}", z, a, T9);
        const auto [bound, data, upperIndex, lowerIndex] = find(z, a, T9);
        if (bound == FRONT || bound == BACK) {
            LOG_TRACE_L2(m_logger, "Derivative is zero for Z={} A={} T9={} (bound: {})", z, a, T9, bound == FRONT ? "FRONT" : "BACK");
            return 0.0; // Derivative is zero at the boundaries
        }
        const auto [T9_high, G_norm_high, T9_low, G_norm_low] = get_interpolation_points(
            upperIndex,
            lowerIndex,
            data.normalized_g_values
        );
        const double slope_g_norm = (G_norm_high - G_norm_low) / (T9_high - T9_low);
        return slope_g_norm * (2.0 * data.ground_state_spin + 1.0);
    }

    bool RauscherThielemannPartitionFunction::supports(
        const int z,
        const int a
    ) const {
        return m_partitionData.contains(make_key(z, a));
    }

    RauscherThielemannPartitionFunction::InterpolationPoints RauscherThielemannPartitionFunction::get_interpolation_points(
        const size_t upper_index,
        const size_t lower_index,
        const std::array<double, 24>& normalized_g_values
    ) {
        const double T_high = RT_TEMPERATURE_GRID_T9[upper_index];
        const double G_norm_high = normalized_g_values[upper_index];
        const double T_low = RT_TEMPERATURE_GRID_T9[lower_index];
        const double G_norm_low = normalized_g_values[lower_index];
        return {T_high, G_norm_high, T_low, G_norm_low};
    }

    RauscherThielemannPartitionFunction::IdentifiedIsotope RauscherThielemannPartitionFunction::find(
        const int z,
        const int a,
        const double T9
    ) const {
        const auto key = make_key(z, a);
        const auto it = m_partitionData.find(key);

        if (it == m_partitionData.end()) {
            LOG_ERROR(m_logger, "Rauscher-Thielemann partition function data for Z={} A={} not found.", z, a);
            throw std::out_of_range("Partition function data not found for Z=" + std::to_string(z) + " A=" + std::to_string(a));
        }

        const IsotopeData& data = it->second;
        const auto upper_it = std::ranges::lower_bound(RT_TEMPERATURE_GRID_T9, T9);
        Bounds bound;
        if (upper_it == RT_TEMPERATURE_GRID_T9.begin()) {
            bound = FRONT; // T9 is below the first grid point
        } else if (upper_it == RT_TEMPERATURE_GRID_T9.end()) {
            bound = BACK; // T9 is above the last grid point
        } else {
            bound = MIDDLE; // T9 is within the grid
        }
        const size_t upper_index = std::distance(RT_TEMPERATURE_GRID_T9.begin(), upper_it);
        const size_t lower_index = upper_index - 1;
        return {bound, data, upper_index, lower_index};
    }

    constexpr int RauscherThielemannPartitionFunction::make_key(
        const int z,
        const int a
    ) {
        return z * 1000 + a; // Simple key generation for Z and A
    }
}
