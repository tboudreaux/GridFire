#pragma once

#include <map>
#include <string>
#include <ostream>
#include <sstream>
#include "fourdst/composition/composition.h"

namespace gridfire::engine {

    /**
     * @enum PrimingReportStatus
     * @brief Enumerates outcome codes for a network priming operation.
     *
     * These status codes indicate the reason for success or failure of the priming process:
     *   - FULL_SUCCESS: Priming completed successfully with all species processed.
     *   - NO_SPECIES_TO_PRIME: There were no species eligible for priming.
     *   - MAX_ITERATIONS_REACHED: The algorithm reached its iteration limit without converging.
     *   - FAILED_TO_FINALIZE_COMPOSITION: Unable to build a valid Composition object at end.
     *   - FAILED_TO_FIND_CREATION_CHANNEL: No reaction path found to create the priming species.
     *   - FAILED_TO_FIND_PRIMING_REACTIONS: No reactions containing the priming species were found.
     *   - BASE_NETWORK_TOO_SHALLOW: The provided base network depth was insufficient for priming.
     *
     * @see PrimingReport for data associated with each status.
     */
    enum class PrimingReportStatus {
        SUCCESS,
        ALREADY_PRIMED,
        SOLVER_FAILURE,
    };

    /**
     * @brief Mapping from PrimingReportStatus codes to human-readable strings.
     *
     * Used when formatting or logging the priming status. No preconditions.
     * The map contains entries for all PrimingReportStatus values.
     */
    inline std::map<PrimingReportStatus, std::string> PrimingReportStatusStrings = {
        {PrimingReportStatus::SUCCESS, "SUCCESS"},
        {PrimingReportStatus::SOLVER_FAILURE, "SOLVER_FAILURE"},
        {PrimingReportStatus::ALREADY_PRIMED, "ALREADY_PRIMED"},
    };

    /**
     * @struct PrimingReport
     * @brief Captures the result of a network priming operation.
     *
     * Contains the finalized composition, the changes in mass fractions for species,
     * a success flag, and a detailed status code.
     *
     * @pre A priming algorithm must be executed before constructing a valid report.
     * @post All fields reflect the outcome of the priming run.
     */
    struct PrimingReport {
        /** Finalized composition after priming. */
        fourdst::composition::Composition primedComposition;
        /** True if priming completed without error. */
        bool success;
        /** Detailed status code indicating the result. */
        PrimingReportStatus status;

        /**
         * @brief Serialize the report to a stream.
         *
         * Formats the success flag and status string into the output stream.
         * @param os Output stream to write to.
         * @param report Report to serialize.
         * @return Reference to the modified output stream.
         * @post The stream contains a textual representation of the report.
         */
        friend std::ostream& operator<<(std::ostream& os, const PrimingReport& report) {
            std::stringstream ss;
            const std::string successStr = report.success ? "true" : "false";
            ss << "PrimingReport(success=" << successStr
               << ", status=" << PrimingReportStatusStrings[report.status] << ")";
            return os << ss.str();
        }
    };

    /**
     * @enum SpeciesStatus
     * @brief Enumerates the status of a species in the simulation.
     *
     * These status codes indicate the current state of a species:
     *   - ACTIVE: The species is actively participating in reactions.
     *   - EQUILIBRIUM: The species is in equilibrium and not changing concentration.
     *   - INACTIVE_FLOW: The species is present but not currently flowing.
     *   - NOT_PRESENT: The species is not present in the system.
     */
    enum class SpeciesStatus {
        ACTIVE,
        EQUILIBRIUM,
        INACTIVE_FLOW,
        NOT_PRESENT
    };

    /**
     * @brief Convert a SpeciesStatus enum value to its string representation.
     *
     * @param status The SpeciesStatus value to convert.
     * @return A string representing the SpeciesStatus.
     */
    inline std::string SpeciesStatus_to_string(const SpeciesStatus status) {
        switch (status) {
            case SpeciesStatus::ACTIVE:
                return "ACTIVE";
            case SpeciesStatus::EQUILIBRIUM:
                return "EQUILIBRIUM";
            case SpeciesStatus::INACTIVE_FLOW:
                return "INACTIVE_FLOW";
            case SpeciesStatus::NOT_PRESENT:
                return "NOT_PRESENT";
            default:
                return "UNKNOWN_STATUS";
        }
    }

}