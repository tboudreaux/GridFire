#pragma once

#include "gridfire/partition/partition_abstract.h"

#include "fourdst/logging/logging.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <utility>
#include <memory>

namespace gridfire::partition {

    /**
     * @class RauscherThielemannPartitionFunction
     * @brief Partition function using Rauscher-Thielemann tabulated normalized G-values.
     *
     * Loads isotope partition data from embedded records and computes values by
     * selecting boundary data or interpolating between grid points on a fixed T9 grid.
     * Implementation in partition_rauscher_thielemann.cpp.
     *
     * @throws std::out_of_range If requested isotope data is missing.
     */
    class RauscherThielemannPartitionFunction final : public PartitionFunction {
    public:
        /**
         * @brief Construct and populate partition data.
         *
         * Reads embedded RT partition data records and fills m_partitionData.
         * @pre Embedded data arrays are available and non-empty.
         * @post m_partitionData contains entries for all isotopes in data.
         */
        RauscherThielemannPartitionFunction();
        /**
         * @brief Evaluate partition function for isotope at temperature.
         *
         * Retrieves boundary or interpolated normalized G-value and scales by (2J+1).
         * @param z Atomic number of the isotope (>=1).
         * @param a Mass number of the isotope (>=z).
         * @param T9 Temperature in units of 10^9 K.
         * @return Dimensionless partition function.
         * @pre supports(z,a) returns true.
         * @post No side effects.
         * @throws std::out_of_range If isotope key not found in m_partitionData.
         */
        double evaluate(int z, int a, double T9) const override;
        /**
         * @brief Evaluate temperature derivative of partition function.
         *
         * Zero at grid extremes; otherwise derivative of linear interpolation.
         * @param z Atomic number (>=1).
         * @param a Mass number (>=z).
         * @param T9 Temperature in 10^9 K.
         * @return d(PartitionFunction)/dT9.
         * @pre supports(z,a) returns true.
         * @post No side effects.
         * @throws std::out_of_range If isotope data is missing.
         */
        double evaluateDerivative(int z, int a, double T9) const override;
        /**
         * @brief Check if partition data exists for given isotope.
         * @param z Atomic number.
         * @param a Mass number.
         * @return true if data available; false otherwise.
         * @post No side effects.
         */
        bool supports(int z, int a) const override;
        /**
         * @brief Get type identifier for this partition function.
         * @return Literal string "RauscherThielemann".
         * @post No side effects.
         */
        std::string type() const override { return "RauscherThielemann"; }
        /**
         * @brief Clone this partition function instance.
         * @return Unique pointer to a copy of this object.
         * @post Caller owns the returned object.
         */
        std::unique_ptr<PartitionFunction> clone() const override {
            return std::make_unique<RauscherThielemannPartitionFunction>(*this);
        }
    private:
        /**
         * @enum Bounds
         * @brief Indicator for temperature grid bound position.
         */
        enum Bounds {
            FRONT,  ///< Below first grid point
            BACK,   ///< Above last grid point
            MIDDLE  ///< Between grid points
        };
    private:
        struct IsotopeData {
            double ground_state_spin;        ///< Spin of the isotope ground state
            std::array<double, 24> normalized_g_values; ///< Normalized G values on RT grid
        };
        struct InterpolationPoints {
            double T9_high;     ///< Upper temperature bound
            double G_norm_high; ///< Normalized G at upper bound
            double T9_low;      ///< Lower temperature bound
            double G_norm_low;  ///< Normalized G at lower bound
        };
        struct IdentifiedIsotope {
            Bounds bound;
            const IsotopeData& data;
            size_t upperIndex;
            size_t lowerIndex;
        };
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::unordered_map<int, IsotopeData> m_partitionData; ///< Map of isotope key to data
    private:
        /**
         * @brief Get interpolation points from normalized G array.
         * @param upper_index Index of upper grid point.
         * @param lower_index Index of lower grid point.
         * @param normalized_g_values Array of normalized G values.
         * @return InterpolationPoints containing bounds and G values.
         */
        static InterpolationPoints get_interpolation_points(
            const size_t upper_index,
            const size_t lower_index,
            const std::array<double, 24>& normalized_g_values
        );
        /**
         * @brief Identify isotope entry and grid indices for given T9.
         * @param z Atomic number of isotope.
         * @param a Mass number of isotope.
         * @param T9 Temperature in 10^9 K.
         * @return IdentifiedIsotope with data reference and indices.
         * @throws std::out_of_range If isotope not found in m_partitionData.
         */
        IdentifiedIsotope find(int z, int a, double T9) const;
        /**
         * @brief Generate integer key for isotope (z,a).
         * @param z Atomic number.
         * @param a Mass number (<1000).
         * @return Key computed as z*1000 + a.
         */
        static constexpr int make_key(int z, int a);
    };

}
