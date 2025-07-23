#pragma once

#include <map>
#include <string>
#include <ranges>

namespace gridfire {
    enum class PrimingReportStatus {
        FULL_SUCCESS = 0,
        NO_SPECIES_TO_PRIME = 1,
        MAX_ITERATIONS_REACHED = 2,
        FAILED_TO_FINALIZE_COMPOSITION = 3,
        FAILED_TO_FIND_CREATION_CHANNEL = 4,
        FAILED_TO_FIND_PRIMING_REACTIONS = 5,
        BASE_NETWORK_TOO_SHALLOW = 6
    };

    inline std::map<PrimingReportStatus, std::string> PrimingReportStatusStrings = {
        {PrimingReportStatus::FULL_SUCCESS, "Full Success"},
        {PrimingReportStatus::NO_SPECIES_TO_PRIME, "No Species to Prime"},
        {PrimingReportStatus::MAX_ITERATIONS_REACHED, "Max Iterations Reached"},
        {PrimingReportStatus::FAILED_TO_FINALIZE_COMPOSITION, "Failed to Finalize Composition"},
        {PrimingReportStatus::FAILED_TO_FIND_CREATION_CHANNEL, "Failed to Find Creation Channel"},
        {PrimingReportStatus::FAILED_TO_FIND_PRIMING_REACTIONS, "Failed to Find Priming Reactions"},
        {PrimingReportStatus::BASE_NETWORK_TOO_SHALLOW, "Base Network Too Shallow"}
    };

    struct PrimingReport {
        fourdst::composition::Composition primedComposition;
        std::vector<std::pair<fourdst::atomic::Species, double>> massFractionChanges; ///< Species and their destruction/creation rates
        bool success;
        PrimingReportStatus status;

        friend std::ostream& operator<<(std::ostream& os, const PrimingReport& report) {
            std::stringstream ss;
            const std::string successStr = report.success ? "true" : "false";
            ss << "PrimingReport(success=" << successStr
               << ", status=" << PrimingReportStatusStrings[report.status] << ")";
            return os << ss.str();
        }
    };
}