#pragma once

#include <variant>

namespace gridfire {

    /**
     * @enum NetworkBuildDepth
     * @brief Specifies supported depths for building the reaction network.
     *
     * Values:
     *   - Full: Build the complete network (infinite depth).
     *   - Shallow: Build only direct reactions (depth = 1).
     *   - SecondOrder: Include reactions up to second order (depth = 2).
     *   - ThirdOrder: Include reactions up to third order (depth = 3).
     *   - FourthOrder: Include reactions up to fourth order (depth = 4).
     *   - FifthOrder: Include reactions up to fifth order (depth = 5).
     *
     * @note For custom build depths, see BuildDepthType.
     */
    enum class NetworkBuildDepth {
        Full = -1,
        Shallow = 1,
        SecondOrder = 2,
        ThirdOrder = 3,
        FourthOrder = 4,
        FifthOrder = 5
    };

    /**
     * @typedef BuildDepthType
     * @brief Variant specifying either a predefined NetworkBuildDepth or a custom integer depth.
     *
     * @pre If using the integer alternative, the value must be >= 0 or -1 to indicate a full build.
     * @post The network builder will interpret and apply the specified depth to control reaction expansion.
     */
    using BuildDepthType = std::variant<NetworkBuildDepth, int>;

}
