#include "gridfire/partition/partition_ground.h"

#include <ranges>

#include "fourdst/logging/logging.h"
#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/species.h"

#include "quill/LogMacros.h"

namespace gridfire::partition {
    GroundStatePartitionFunction::GroundStatePartitionFunction() {
        for (const auto &isotope: fourdst::atomic::species | std::views::values) {
            m_ground_state_spin[make_key(isotope.z(), isotope.a())] = isotope.spin();
        }
    }

    double GroundStatePartitionFunction::evaluate(
        const int z,
        const int a,
        const double T9
    ) const {
        LOG_TRACE_L2(m_logger, "Evaluating ground state partition function for Z={} A={} T9={}", z, a, T9);
        const int key = make_key(z, a);
        const double spin = m_ground_state_spin.at(key);
        return (2.0 * spin) + 1.0;
    }

    double GroundStatePartitionFunction::evaluateDerivative(
        const int z,
        const int a,
        const double T9
    ) const {
        LOG_TRACE_L2(m_logger, "Evaluating derivative of ground state partition function for Z={} A={} T9={}", z, a, T9);
        return 0.0;
    }

    bool GroundStatePartitionFunction::supports(
        const int z,
        const int a
    ) const {
        return m_ground_state_spin.contains(make_key(z, a));
    }
    constexpr int GroundStatePartitionFunction::make_key(
        const int z,
        const int a
    ) {
        return z * 1000 + a; // Simple key generation for Z and A
    }
}
