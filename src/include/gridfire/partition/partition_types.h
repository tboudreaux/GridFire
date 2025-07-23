#pragma once

#include <unordered_map>
#include <string>

namespace gridfire::partition {
    enum BasePartitionType {
        RauscherThielemann, ///< Rauscher-Thielemann partition function
        GroundState,        ///< Ground state partition function
    };

    inline std::unordered_map<BasePartitionType, std::string> basePartitionTypeToString = {
        {RauscherThielemann, "RauscherThielemann"},
        {GroundState, "GroundState"}
    };

    inline std::unordered_map<std::string, BasePartitionType> stringToBasePartitionType = {
        {"RauscherThielemann", RauscherThielemann},
        {"GroundState", GroundState}
    };
}