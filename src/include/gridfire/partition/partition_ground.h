#pragma once

#include "gridfire/partition/partition_abstract.h"

#include "fourdst/logging/logging.h"

#include <unordered_map>
#include <memory>

#include "quill/Logger.h"

namespace gridfire::partition {
    /**
     * @class GroundStatePartitionFunction
     * @brief Partition function implementation for nuclear ground states.
     *
     * Computes the partition function as (2J + 1) based on the ground state spin J of each isotope.
     * The temperature derivative is always zero. Ground state spins are loaded from the
     * fourdst::atomic::species registry at construction.
     * @see partition_ground.cpp for implementation details.
     */
     class GroundStatePartitionFunction final : public PartitionFunction {
     public:
        /**
         * @brief Construct and populate the ground state spin map.
+         *
         * Loads spins for all isotopes from the atomic species registry into m_ground_state_spin.
         * @pre atomic::species registry is initialized and non-empty.
         * @post m_ground_state_spin contains entries for each isotope.
+         */
         GroundStatePartitionFunction();
        /**
         * @brief Evaluate the ground state partition function.
+         *
         * @param z Proton number (atomic number) of the isotope; must be >= 1.
         * @param a Mass number of the isotope; must be >= z.
         * @param T9 Temperature in units of 10^9 K; unused for ground state.
         * @pre supports(z,a) returns true.
         * @post No side effects.
         * @return Dimensionless partition function value = 2*spin + 1.
         * @throws std::out_of_range If the isotope key is not found in m_ground_state_spin.
+         */
         double evaluate(
             const int z,
             const int a,
             const double T9
         ) const override;
        /**
         * @brief Evaluate the temperature derivative of the ground state partition function.
+         *
         * Always returns zero as ground state has no temperature dependence.
         * @param z Proton number of the isotope; must be supported.
         * @param a Mass number of the isotope; must be supported.
         * @param T9 Temperature in units of 10^9 K; unused.
         * @pre supports(z,a) returns true.
         * @post No side effects.
         * @return Zero.
         * @throws std::out_of_range If the isotope key is not found.
+         */
         double evaluateDerivative(
             const int z,
             const int a,
             const double T9
         ) const override;
        /**
         * @brief Check if ground state data exists for the given isotope.
+         *
         * @param z Proton number of the isotope.
         * @param a Mass number of the isotope.
         * @return True if m_ground_state_spin contains the key; false otherwise.
         * @post No side effects.
+         */
         bool supports(
             const int z,
             const int a
         ) const override;
        /**
         * @brief Get the type identifier of this partition function.
         * @return The string literal "GroundState".
         * @post No side effects.
+         */
         std::string type() const override { return "GroundState"; }
        /**
         * @brief Create a deep copy of this partition function.
         * @return Unique_ptr to a new GroundStatePartitionFunction cloned from this object.
         * @post Caller owns the returned instance.
+         */
         std::unique_ptr<PartitionFunction> clone() const override {
             return std::make_unique<GroundStatePartitionFunction>(*this);
         }
     private:
         quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
         std::unordered_map<int, double> m_ground_state_spin;
        /**
         * @brief Generate a unique lookup key for an isotope.
+         *
         * Combines atomic number z and mass number a into a single integer.
         * @param z Proton number of the isotope.
         * @param a Mass number of the isotope; should be < 1000 to avoid collisions.
         * @pre a < 1000.
         * @return Integer key = z * 1000 + a.
+         */
         static constexpr int make_key(
             const int z,
         const int a);
     };

 }
