#pragma once

#include <map>
#include <string>
#include <ranges>
// Required for PrimingReport fields and streaming
#include <vector>
#include <utility>
#include <ostream>
#include <sstream>
#include "fourdst/composition/composition.h"
#include "fourdst/composition/atomicSpecies.h"

namespace gridfire {

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
        FULL_SUCCESS = 0,
        NO_SPECIES_TO_PRIME = 1,
        MAX_ITERATIONS_REACHED = 2,
        FAILED_TO_FINALIZE_COMPOSITION = 3,
        FAILED_TO_FIND_CREATION_CHANNEL = 4,
        FAILED_TO_FIND_PRIMING_REACTIONS = 5,
        BASE_NETWORK_TOO_SHALLOW = 6
    };

    /**
     * @brief Mapping from PrimingReportStatus codes to human-readable strings.
     *
     * Used when formatting or logging the priming status. No preconditions.
     * The map contains entries for all PrimingReportStatus values.
     */
    inline std::map<PrimingReportStatus, std::string> PrimingReportStatusStrings = {
        {PrimingReportStatus::FULL_SUCCESS, "Full Success"},
        {PrimingReportStatus::NO_SPECIES_TO_PRIME, "No Species to Prime"},
        {PrimingReportStatus::MAX_ITERATIONS_REACHED, "Max Iterations Reached"},
        {PrimingReportStatus::FAILED_TO_FINALIZE_COMPOSITION, "Failed to Finalize Composition"},
        {PrimingReportStatus::FAILED_TO_FIND_CREATION_CHANNEL, "Failed to Find Creation Channel"},
        {PrimingReportStatus::FAILED_TO_FIND_PRIMING_REACTIONS, "Failed to Find Priming Reactions"},
        {PrimingReportStatus::BASE_NETWORK_TOO_SHALLOW, "Base Network Too Shallow"}
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
        /**
         * List of pairs (species, rate change) representing destruction (<0)
         * or creation (>0) rates of species during priming.
         */
        std::vector<std::pair<fourdst::atomic::Species, double>> massFractionChanges;
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

}