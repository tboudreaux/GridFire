#pragma once

#include <variant>

namespace gridfire {
    enum class NetworkBuildDepth {
        Full = -1,
        Shallow = 1,
        SecondOrder = 2,
        ThirdOrder = 3,
        FourthOrder = 4,
        FifthOrder = 5
    };

    using BuildDepthType = std::variant<NetworkBuildDepth, int>;
}
