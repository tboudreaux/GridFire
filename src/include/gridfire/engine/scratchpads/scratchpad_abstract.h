/**
 * @file scratchpad_abstract.h
 * @brief Abstract base class for scratchpad memory used during engine computations.
 *
 * This header defines the AbstractScratchPad interface, which provides a common
 * contract for temporary working memory (scratchpads) used by computational engines.
 * Scratchpads are designed to store intermediate results, cached computations, or
 * pre-allocated buffers that can be reused across multiple computational steps,
 * improving performance by avoiding repeated memory allocations.
 *
 * @par Purpose
 * The scratchpad pattern allows engines to:
 * - Pre-allocate working memory once and reuse it across iterations
 * - Store intermediate computational results between solver steps
 * - Enable efficient cloning for parallel execution contexts
 * - Provide a type-erased interface for heterogeneous scratchpad management
 *
 * @par Examples
 * @code{.cpp}
 * // Define a concrete scratchpad for a specific engine
 * class MyScratchPad : public gridfire::engine::scratch::AbstractScratchPad {
 * public:
 *     MyScratchPad() : initialized_(false) {}
 *
 *     void initialize(size_t size) {
 *         buffer_.resize(size);
 *         initialized_ = true;
 *     }
 *
 *     [[nodiscard]] bool is_initialized() const override {
 *         return initialized_;
 *     }
 *
 *     [[nodiscard]] std::unique_ptr<AbstractScratchPad> clone() const override {
 *         auto copy = std::make_unique<MyScratchPad>();
 *         copy->buffer_ = buffer_;
 *         copy->initialized_ = initialized_;
 *         return copy;
 *     }
 *
 * private:
 *     std::vector<double> buffer_;
 *     bool initialized_;
 * };
 *
 * // Usage in an engine context
 * auto scratch = std::make_unique<MyScratchPad>();
 * scratch->initialize(1024);
 *
 * if (scratch->is_initialized()) {
 *     // Use scratchpad for computations
 *     auto parallel_scratch = scratch->clone();  // Clone for parallel worker
 * }
 * @endcode
 */

#pragma once
#include <memory>

namespace gridfire::engine::scratch {

/**
 * @brief Abstract base struct for engine scratchpad memory.
 *
 * AbstractScratchPad defines the interface for temporary working memory
 * containers used by computational engines. Implementations should provide
 * storage for intermediate results, cached values, or pre-allocated buffers
 * that persist across multiple computational steps.
 *
 * This interface enables polymorphic handling of different scratchpad types
 * while ensuring proper resource management through virtual destruction and
 * deep cloning capabilities.
 *
 * @par Thread Safety
 * This interface is **not thread-safe** by design. Scratchpads are intended
 * to be used as thread-local working memory. Each thread should operate on
 * its own independent scratchpad instance. Use the clone() method to create
 * separate copies for each thread in parallel execution contexts. Sharing a
 * single scratchpad instance across multiple threads without external
 * synchronization will result in undefined behavior.
 */
struct AbstractScratchPad {
    /**
     * @brief Virtual destructor for proper cleanup of derived classes.
     *
     * Ensures that resources held by concrete scratchpad implementations
     * are properly released when the scratchpad is destroyed through a
     * base class pointer.
     */
    virtual ~AbstractScratchPad() = default;

    /**
     * @brief Check whether the scratchpad has been properly initialized.
     *
     * Derived classes should return true only after all necessary memory
     * allocations and setup operations have been completed successfully.
     *
     * @return true if the scratchpad is initialized and ready for use.
     * @return false if the scratchpad has not been initialized or initialization failed.
     *
     * @par Examples
     * @code{.cpp}
     * auto scratch = create_scratchpad();
     * if (!scratch->is_initialized()) {
     *     throw std::runtime_error("Scratchpad not ready for computation");
     * }
     * @endcode
     */
    [[nodiscard]] virtual bool is_initialized() const = 0;

    /**
     * @brief Create a deep copy of this scratchpad.
     *
     * Produces an independent clone of the scratchpad, including all internal
     * state and allocated memory. This is essential for parallel execution
     * scenarios where each thread requires its own working memory.
     *
     * @return A unique pointer to a newly allocated copy of this scratchpad.
     *
     * @note The returned clone should be fully independent; modifications to
     *       the clone must not affect the original, and vice versa.
     *
     * @par Examples
     * @code{.cpp}
     * std::unique_ptr<AbstractScratchPad> original = create_scratchpad();
     *
     * // Create independent copies for parallel workers
     * std::vector<std::unique_ptr<AbstractScratchPad>> worker_scratches;
     * for (int i = 0; i < num_threads; ++i) {
     *     worker_scratches.push_back(original->clone());
     * }
     * @endcode
     */
    [[nodiscard]] virtual std::unique_ptr<AbstractScratchPad> clone() const = 0;
};

} // namespace gridfire::engine::scratch
