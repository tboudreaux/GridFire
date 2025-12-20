/**
 * @file blob.h
 * @brief Container class for managing multiple scratchpad instances.
 *
 * This header defines the StateBlob class, which serves as a centralized
 * registry for managing multiple scratchpad instances used by computational
 * engines. It provides type-safe enrollment, retrieval, and cloning of
 * scratchpads using compile-time type checking via C++20 concepts.
 *
 * @par Purpose
 * The StateBlob provides:
 * - A fixed-size array of scratchpad slots indexed by ScratchPadType
 * - Type-safe enrollment ensuring one instance per scratchpad type
 * - Compile-time verified retrieval with optional initialization checks
 * - Deep cloning of all enrolled scratchpads for parallel execution
 * - Status tracking for each scratchpad slot
 *
 * @par Examples
 * @code{.cpp}
 * #include "gridfire/engine/scratchpads/blob.h"
 * #include "gridfire/engine/scratchpads/engine_graph_scratchpad.h"
 *
 * using namespace gridfire::engine::scratch;
 *
 * // Create a StateBlob and enroll scratchpads
 * StateBlob blob;
 * blob.enroll<GraphEngineScratchPad>();
 *
 * // Retrieve a scratchpad (returns std::expected)
 * auto result = blob.get<GraphEngineScratchPad>();
 * if (result.has_value()) {
 *     GraphEngineScratchPad* scratch = result.value();
 *     scratch->initialize(engine);
 * }
 *
 * // Retrieve with initialization check
 * auto checked = blob.get<GraphEngineScratchPad, true>();
 * if (!checked.has_value()) {
 *     if (checked.error() == StateBlob::Error::SCRATCHPAD_NOT_INITIALIZED) {
 *         // Handle uninitialized scratchpad
 *     }
 * }
 *
 * // Clone for parallel execution
 * auto worker_blob = blob.clone_structure();
 * @endcode
 *
 * @par Thread Safety
 * The StateBlob class is **not thread-safe**. Each thread should have its own
 * StateBlob instance. Use clone_structure() to create independent copies for
 * parallel workers. The cloned blob contains deep copies of all enrolled
 * scratchpads.
 *
 * @see AbstractScratchPad
 * @see ScratchPadType
 */

#pragma once

#include "gridfire/engine/scratchpads/scratchpad_abstract.h"
#include "gridfire/engine/scratchpads/types.h"
#include "gridfire/exceptions/error_scratchpad.h"

#include <unordered_map>
#include <memory>
#include <expected>
#include <unordered_set>

namespace gridfire::engine::scratch {

/**
 * @brief Concept that constrains types to valid scratchpad implementations.
 *
 * A type satisfies IsScratchPad if:
 * - It derives from AbstractScratchPad
 * - It has a static ID member convertible to ScratchPadType
 *
 * @tparam T The type to check against the concept.
 *
 * @par Examples
 * @code{.cpp}
 * // This will compile only if MyScratchPad satisfies IsScratchPad
 * template <IsScratchPad T>
 * void process_scratchpad(T& scratch) {
 *     // Use scratch...
 * }
 * @endcode
 */
template <typename T>
concept IsScratchPad = std::is_base_of_v<AbstractScratchPad, T>
                       && requires { { T::ID } -> std::convertible_to<ScratchPadType>; };


/**
 * @brief Container for managing a collection of typed scratchpad instances.
 *
 * StateBlob provides a centralized registry for scratchpads used by engines.
 * It uses a fixed-size array indexed by ScratchPadType for O(1) access, with
 * compile-time type safety enforced through the IsScratchPad concept.
 *
 * The blob supports:
 * - Enrolling new scratchpad types (one instance per type)
 * - Type-safe retrieval with optional initialization validation
 * - Runtime retrieval by ScratchPadType enum value
 * - Deep cloning for parallel execution contexts
 * - Status queries for monitoring scratchpad states
 *
 * @par Thread Safety
 * This class is **not thread-safe**. Each thread should have its own StateBlob
 * instance. Use clone_structure() to create independent copies for parallel
 * workers. Concurrent access to the same StateBlob instance requires external
 * synchronization.
 */
class StateBlob {
public:
    /**
     * @brief Error codes for scratchpad operations.
     */
    enum class Error : uint8_t {
        SCRATCHPAD_NOT_FOUND,       ///< Requested scratchpad type is not enrolled.
        SCRATCHPAD_BAD_CAST,        ///< Dynamic cast to requested type failed.
        SCRATCHPAD_NOT_INITIALIZED, ///< Scratchpad exists but is not initialized.
        SCRATCHPAD_TYPE_COLLISION,  ///< Attempted to enroll duplicate type.
        SCRATCHPAD_OUT_OF_BOUNDS,   ///< ScratchPadType index exceeds array bounds.
        SCRATCHPAD_UNKNOWN_ERROR    ///< Unspecified error condition.
    };

    /**
     * @brief Convert an Error enum value to a human-readable string.
     *
     * @param error The error code to convert.
     * @return A string representation of the error.
     *
     * @par Examples
     * @code{.cpp}
     * auto result = blob.get<MyScratchPad>();
     * if (!result.has_value()) {
     *     std::cerr << "Error: " << StateBlob::error_to_string(result.error());
     * }
     * @endcode
     */
    static std::string error_to_string(const Error error) {
        switch (error) {
            case Error::SCRATCHPAD_NOT_FOUND:
                return "SCRATCHPAD_NOT_FOUND";
            case Error::SCRATCHPAD_BAD_CAST:
                return "SCRATCHPAD_BAD_CAST";
            case Error::SCRATCHPAD_NOT_INITIALIZED:
                return "SCRATCHPAD_NOT_INITIALIZED";
            case Error::SCRATCHPAD_TYPE_COLLISION:
                return "SCRATCHPAD_TYPE_COLLISION";
            case Error::SCRATCHPAD_OUT_OF_BOUNDS:
                return "SCRATCHPAD_OUT_OF_BOUNDS";
            default:
                return "SCRATCHPAD_UNKNOWN_ERROR";
        }
    }

    /**
     * @brief Status codes for scratchpad slots.
     */
    enum class ScratchPadStatus : uint8_t {
        NOT_ENROLLED,            ///< No scratchpad has been enrolled for this slot.
        ENROLLED_NOT_INITIALIZED,///< Scratchpad enrolled but not yet initialized.
        ENROLLED_INITIALIZED     ///< Scratchpad enrolled and fully initialized.
    };

public:
    /// @brief Default constructor.
    StateBlob() = default;

    /// @brief Default destructor.
    ~StateBlob() = default;

    /**
     * @brief Enroll a new scratchpad type into the blob.
     *
     * Creates a new instance of the specified scratchpad type and registers it
     * in the appropriate slot. Only one instance per type is allowed.
     *
     * @tparam CTX The scratchpad type to enroll (must satisfy IsScratchPad).
     *
     * @throws exceptions::ScratchPadError if a scratchpad of this type is already enrolled.
     *
     * @par Examples
     * @code{.cpp}
     * StateBlob blob;
     * blob.enroll<GraphEngineScratchPad>();
     * blob.enroll<AdaptiveEngineViewScratchPad>();
     *
     * // This would throw - duplicate enrollment
     * // blob.enroll<GraphEngineScratchPad>();
     * @endcode
     */
    template <IsScratchPad CTX>
    void enroll() {
        constexpr auto index = static_cast<size_t>(CTX::ID);
        static_assert(index < MAX_SCRATCHPADS, "ScratchPadType ID exceeds (maximum) allowed scratchpads.");
        if (scratchpad_enrolled_flags[index]) {
            throw exceptions::ScratchPadError("ScratchPad of this type has already been enrolled. Only one instance per type is allowed.");
        }
        scratchpads[index] = std::make_unique<CTX>();
        scratchpad_enrolled_flags[index] = true;
    }

    /**
     * @brief Retrieve a scratchpad by type.
     *
     * Returns a pointer to the enrolled scratchpad of the specified type.
     * In debug builds, performs a dynamic_cast for type safety; in release
     * builds, uses static_cast for performance.
     *
     * @tparam CTX The scratchpad type to retrieve (must satisfy IsScratchPad).
     *
     * @return std::expected containing the scratchpad pointer, or an Error if not found/invalid.
     *
     * @par Examples
     * @code{.cpp}
     * auto result = blob.get<GraphEngineScratchPad>();
     * if (result.has_value()) {
     *     GraphEngineScratchPad* scratch = result.value();
     *     // Use scratch...
     * } else {
     *     // Handle error
     *     std::cerr << StateBlob::error_to_string(result.error());
     * }
     * @endcode
     */
    template <IsScratchPad CTX>
    std::expected<CTX*, Error> get() const {
        constexpr auto index = static_cast<size_t>(CTX::ID);
        static_assert(index < MAX_SCRATCHPADS, "ScratchPadType ID exceeds maximum allowed scratchpads.");
        AbstractScratchPad* scratchpad = scratchpads[index].get();
        if (!scratchpad) {
            return std::unexpected<Error>(Error::SCRATCHPAD_NOT_FOUND);
        }
#if !defined(NDEBUG)
        if (auto* cast_ptr = dynamic_cast<CTX*>(scratchpad)) {
            return cast_ptr;
        } else {
            return std::unexpected<Error>(Error::SCRATCHPAD_BAD_CAST);
        }
#else
        return static_cast<CTX*>(scratchpad);
#endif
    }

    /**
     * @brief Retrieve a scratchpad by type with optional initialization check.
     *
     * Returns a pointer to the enrolled scratchpad of the specified type.
     * When MUST_BE_INITIALIZED is true, also verifies that the scratchpad
     * has been initialized before returning it.
     *
     * @tparam CTX The scratchpad type to retrieve (must satisfy IsScratchPad).
     * @tparam MUST_BE_INITIALIZED If true, returns an error for uninitialized scratchpads.
     *
     * @return std::expected containing the scratchpad pointer, or an Error if not found/invalid/uninitialized.
     *
     * @par Examples
     * @code{.cpp}
     * // Get only if initialized
     * auto result = blob.get<GraphEngineScratchPad, true>();
     * if (!result.has_value()) {
     *     if (result.error() == StateBlob::Error::SCRATCHPAD_NOT_INITIALIZED) {
     *         // Need to initialize first
     *     }
     * }
     * @endcode
     */
    template <IsScratchPad CTX, bool MUST_BE_INITIALIZED>
    std::expected<CTX*, Error> get() const {
        constexpr auto index = static_cast<size_t>(CTX::ID);
        static_assert(index < MAX_SCRATCHPADS, "ScratchPadType ID exceeds maximum allowed scratchpads.");
        AbstractScratchPad* scratchpad = scratchpads[index].get();
        if (!scratchpad) {
            return std::unexpected<Error>(Error::SCRATCHPAD_NOT_FOUND);
        }
#if !defined(NDEBUG)
        if (auto* cast_ptr = dynamic_cast<CTX*>(scratchpad)) {
            if constexpr (MUST_BE_INITIALIZED) {
                if (!cast_ptr->is_initialized()) {
                    return std::unexpected<Error>(Error::SCRATCHPAD_NOT_INITIALIZED);
                }
            }
            return cast_ptr;
        } else {
            return std::unexpected<Error>(Error::SCRATCHPAD_BAD_CAST);
        }
#else
        CTX* cast_ptr = static_cast<CTX*>(scratchpad);
        if constexpr (MUST_BE_INITIALIZED) {
            if (!cast_ptr->is_initialized()) {
                return std::unexpected<Error>(Error::SCRATCHPAD_NOT_INITIALIZED);
            }
        }
        return cast_ptr;
#endif
    }

    /**
     * @brief Retrieve a scratchpad by runtime ScratchPadType value.
     *
     * Returns a pointer to the abstract base class for the scratchpad at
     * the specified type index. Useful when the concrete type is not known
     * at compile time.
     *
     * @param type The ScratchPadType enum value identifying the scratchpad.
     *
     * @return std::expected containing the AbstractScratchPad pointer, or an Error.
     *
     * @par Examples
     * @code{.cpp}
     * ScratchPadType type = ScratchPadType::GRAPH_ENGINE_SCRATCHPAD;
     * auto result = blob.get(type);
     * if (result.has_value()) {
     *     AbstractScratchPad* scratch = result.value();
     *     if (scratch->is_initialized()) {
     *         // Use scratch...
     *     }
     * }
     * @endcode
     */
    [[nodiscard]] std::expected<AbstractScratchPad*, Error> get(const ScratchPadType type) const { // NOLINT(*-convert-member-functions-to-static)
        const auto index = static_cast<size_t>(type);
        if (index >= MAX_SCRATCHPADS) {
            return std::unexpected<Error>(Error::SCRATCHPAD_OUT_OF_BOUNDS);
        }
        AbstractScratchPad* scratchpad = scratchpads[index].get();
        if (!scratchpad) {
            return std::unexpected<Error>(Error::SCRATCHPAD_NOT_FOUND);
        }
        return scratchpad;
    }

    /**
     * @brief Create a deep copy of this blob with all enrolled scratchpads.
     *
     * Clones the blob structure and all enrolled scratchpads using their
     * clone() methods. The resulting blob is independent and can be used
     * in a separate thread.
     *
     * @return A unique pointer to the cloned StateBlob.
     *
     * @par Examples
     * @code{.cpp}
     * StateBlob blob;
     * blob.enroll<GraphEngineScratchPad>();
     *
     * // Initialize the scratchpad
     * auto scratch = blob.get<GraphEngineScratchPad>().value();
     * scratch->initialize(engine);
     *
     * // Clone for parallel workers
     * std::vector<std::unique_ptr<StateBlob>> worker_blobs;
     * for (int i = 0; i < num_threads; ++i) {
     *     worker_blobs.push_back(blob.clone_structure());
     * }
     * @endcode
     */
    [[nodiscard]] std::unique_ptr<StateBlob> clone_structure() const { // NOLINT(*-convert-member-functions-to-static)
        auto new_blob = std::make_unique<StateBlob>();
        for (size_t i = 0; i < MAX_SCRATCHPADS; ++i) {
            if (scratchpad_enrolled_flags[i]) {
                new_blob->scratchpads[i] = scratchpads[i]->clone();
                new_blob->scratchpad_enrolled_flags[i] = true;
            }
        }
        return new_blob;
    }

    /**
     * @brief Get the set of all registered scratchpad types.
     *
     * @return An unordered set of ScratchPadType values for all enrolled scratchpads.
     *
     * @par Examples
     * @code{.cpp}
     * auto registered = blob.get_registered_scratchpads();
     * for (auto type : registered) {
     *     std::cout << "Registered: " << static_cast<int>(type) << "\n";
     * }
     * @endcode
     */
    [[nodiscard]] std::unordered_set<ScratchPadType> get_registered_scratchpads() const { // NOLINT(*-convert-member-functions-to-static)
        std::unordered_set<ScratchPadType> sset;
        for (size_t i = 0; i < MAX_SCRATCHPADS; ++i) {
            if (scratchpad_enrolled_flags[i]) {
                sset.insert(static_cast<ScratchPadType>(i));
            }
        }
        return sset;
    }

    /**
     * @brief Check if a specific scratchpad type is initialized.
     *
     * @tparam CTX The scratchpad type to check (must satisfy IsScratchPad).
     *
     * @return true if the scratchpad is enrolled and initialized.
     *
     * @throws exceptions::ScratchPadError if the scratchpad type is not enrolled.
     *
     * @par Examples
     * @code{.cpp}
     * if (blob.initialized<GraphEngineScratchPad>()) {
     *     // Safe to use the scratchpad
     * }
     * @endcode
     */
    template <IsScratchPad CTX>
    [[nodiscard]] bool initialized() const {
        constexpr auto index = static_cast<size_t>(CTX::ID);
        static_assert(index < MAX_SCRATCHPADS, "ScratchPadType ID exceeds maximum allowed scratchpads.");
        if (!scratchpad_enrolled_flags[index]) {
            throw exceptions::ScratchPadError("Cannot check initialization status: ScratchPad of this type is not enrolled.");
        }
        return scratchpads[index]->is_initialized();
    }

    /**
     * @brief Get the status of a specific scratchpad type.
     *
     * @tparam CTX The scratchpad type to query (must satisfy IsScratchPad).
     *
     * @return The ScratchPadStatus indicating enrollment and initialization state.
     *
     * @par Examples
     * @code{.cpp}
     * auto status = blob.get_status<GraphEngineScratchPad>();
     * switch (status) {
     *     case StateBlob::ScratchPadStatus::NOT_ENROLLED:
     *         blob.enroll<GraphEngineScratchPad>();
     *         break;
     *     case StateBlob::ScratchPadStatus::ENROLLED_NOT_INITIALIZED:
     *         // Need to initialize
     *         break;
     *     case StateBlob::ScratchPadStatus::ENROLLED_INITIALIZED:
     *         // Ready to use
     *         break;
     * }
     * @endcode
     */
    template <IsScratchPad CTX>
    [[nodiscard]] ScratchPadStatus get_status() const {
        constexpr auto index = static_cast<size_t>(CTX::ID);
        static_assert(index < MAX_SCRATCHPADS, "ScratchPadType ID exceeds maximum allowed scratchpads.");
        if (!scratchpad_enrolled_flags[index]) {
            return ScratchPadStatus::NOT_ENROLLED;
        }
        if (scratchpads[index]->is_initialized()) {
            return ScratchPadStatus::ENROLLED_INITIALIZED;
        } else {
            return ScratchPadStatus::ENROLLED_NOT_INITIALIZED;
        }
    }

    /**
     * @brief Get a map of all scratchpad types to their current status.
     *
     * @return An unordered map from ScratchPadType to ScratchPadStatus for all slots.
     *
     * @par Examples
     * @code{.cpp}
     * auto status_map = blob.get_status_map();
     * for (const auto& [type, status] : status_map) {
     *     if (status == StateBlob::ScratchPadStatus::ENROLLED_NOT_INITIALIZED) {
     *         std::cout << "Scratchpad " << static_cast<int>(type) << " needs initialization\n";
     *     }
     * }
     * @endcode
     */
    [[nodiscard]] std::unordered_map<ScratchPadType, ScratchPadStatus> get_status_map() const { // NOLINT(*-convert-member-functions-to-static)
        std::unordered_map<ScratchPadType, ScratchPadStatus> status_map;
        for (size_t i = 0; i < MAX_SCRATCHPADS; ++i) {
            auto type = static_cast<ScratchPadType>(i);
            if (!scratchpad_enrolled_flags[i]) {
                status_map[type] = ScratchPadStatus::NOT_ENROLLED;
            } else if (scratchpads[i]->is_initialized()) {
                status_map[type] = ScratchPadStatus::ENROLLED_INITIALIZED;
            } else {
                status_map[type] = ScratchPadStatus::ENROLLED_NOT_INITIALIZED;
            }
        }
        return status_map;
    }

private:
    /// @brief Maximum number of scratchpad slots, derived from ScratchPadType::_COUNT.
    static constexpr size_t MAX_SCRATCHPADS = static_cast<size_t>(ScratchPadType::_COUNT);

    /// @brief Array of scratchpad instances indexed by ScratchPadType.
    std::array<std::unique_ptr<AbstractScratchPad>, MAX_SCRATCHPADS> scratchpads{};

    /// @brief Flags indicating which scratchpad slots have been enrolled.
    std::array<bool, MAX_SCRATCHPADS> scratchpad_enrolled_flags{false};
};

} // namespace gridfire::engine::scratch
