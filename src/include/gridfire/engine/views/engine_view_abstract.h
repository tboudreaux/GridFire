#pragma once

#include "gridfire/engine/engine_abstract.h"

/**
 * @file engine_view_abstract.h
 * @brief Abstract interfaces for engine "views" in GridFire.
 *
 * This header defines the abstract base classes and concepts for "views" of reaction network engines.
 * The primary purpose of an EngineView is to enable dynamic or adaptive network topologies
 * (such as species/reaction culling, masking, or remapping) without modifying the underlying
 * physics engine or its implementation. Engine views act as a flexible interface layer,
 * allowing the network structure to be changed at runtime while preserving the core
 * physics and solver logic in the base engine.
 *
 * Typical use cases include:
 *   - Adaptive or reduced networks for computational efficiency.
 *   - Dynamic masking or culling of species/reactions based on runtime criteria.
 *   - Providing a filtered or remapped view of the network for integration or analysis.
 *   - Supporting dynamic network reconfiguration in multi-zone or multi-physics contexts.
 *
 * The base engine types referenced here are defined in @ref engine_abstract.h.
 *
 * @author
 * Emily M. Boudreaux
 */

namespace gridfire {

    /**
     * @brief Concept for types allowed as engine bases in EngineView.
     *
     * This concept restricts template parameters to types derived from either
     * gridfire::Engine or gridfire::DynamicEngine, as defined in engine_abstract.h.
     *
     * Example usage:
     * @code
     * static_assert(EngineType<MyEngine>);
     * @endcode
     */
    template<typename EngineT>
    concept EngineType = std::is_base_of_v<Engine, EngineT> || std::is_base_of_v<DynamicEngine, EngineT>;

    /**
     * @brief Abstract base class for a "view" of a reaction network engine.
     *
     * @tparam EngineT The engine type being viewed (must satisfy EngineType).
     *
     * EngineView provides an interface for accessing an underlying engine instance,
     * while presenting a potentially modified or reduced network structure to the user.
     * This enables dynamic or adaptive network topologies (e.g., culling, masking, or
     * remapping of species and reactions) without altering the core physics engine.
     *
     * Intended usage: Derive from this class to implement a custom view or wrapper
     * that manages a dynamic or adaptive network structure, delegating core calculations
     * to the base engine. The contract is that getBaseEngine() must return a reference
     * to the underlying engine instance, which remains responsible for the full physics.
     *
     * Example (see also AdaptiveEngineView):
     * @code
     * class MyAdaptiveView : public gridfire::EngineView<DynamicEngine> {
     * public:
     *     MyAdaptiveView(DynamicEngine& engine) : engine_(engine) {}
     *     const DynamicEngine& getBaseEngine() const override { return engine_; }
     *     // Implement dynamic masking/culling logic...
     * private:
     *     DynamicEngine& engine_;
     * };
     * @endcode
     *
     * @see gridfire::AdaptiveEngineView for a concrete example of dynamic culling.
     */
    template<EngineType EngineT>
    class EngineView {
    public:
        /**
         * @brief Virtual destructor.
         */
        virtual ~EngineView() = default;

        /**
         * @brief Access the underlying engine instance.
         *
         * @return Const reference to the underlying engine.
         *
         * This method must be implemented by derived classes to provide access
         * to the base engine. The returned reference should remain valid for the
         * lifetime of the EngineView.
         *
         * Example:
         * @code
         * const DynamicEngine& engine = myView.getBaseEngine();
         * @endcode
         */
        virtual const EngineT& getBaseEngine() const = 0;
    };

}