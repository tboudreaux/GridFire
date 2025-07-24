/**
 * @file engine.h
 * @brief Core header for the GridFire reaction network engine module.
 *
 * This module defines the core interfaces and classes for reaction network
 * engines in GridFire. It provides abstract base classes for engines,
 * dynamic engines, and engine views, as well as concrete engine
 * implementations and view implementations.
 *
 * The engine module is designed to support a wide range of reaction network
 * simulations, from simple single-zone calculations to complex multi-zone
 * simulations with adaptive network topologies.
 *
 * @section EngineDesign Engine Design
 *
 * The engine module is built around the following key concepts:
 *
 *  - **Engine:** The base class for all reaction network engines. It defines
 *    the minimal interface for evaluating the right-hand side (dY/dt) and
 *    energy generation rate for a given set of abundances, temperature, and
 *    density.
 *
 *  - **DynamicEngine:** An extension of the Engine class that supports
 *    Jacobian and stoichiometry operations, as well as the ability to
 *    dynamically modify the reaction network.
 *
 *  - **EngineView:** An abstract base class for "views" of reaction network
 *    engines. Engine views provide a way to dynamically or adaptively
 *    modify the network topology without modifying the underlying physics
 *    engine.
 *
 * @section EngineComposition Engine Composition
 *
 * Engines and engine views can be composed to create complex reaction network
 * simulations. For example, an AdaptiveEngineView can be used to dynamically
 * cull species and reactions from a GraphEngine, reducing the computational
 * cost of the simulation.
 *
 * The order in which engines and engine views are composed is important. The
 * base engine should always be the innermost engine, and the engine views
 * should be layered on top of the base engine.
 *
 * @section AvailableEngines Available Engines
 *
 * The engine module provides the following concrete engine implementations:
 *
 *  - **GraphEngine:** A reaction network engine that uses a graph-based
 *    representation of the reaction network. It uses sparse matrices for
 *    efficient storage and computation of the stoichiometry and Jacobian
 *    matrices.
 *
 * @section AvailableViews Available Views
 *
 * The engine module provides the following engine view implementations:
 *
 *  - **AdaptiveEngineView:** An engine view that dynamically adapts the
 *    reaction network based on runtime conditions. It culls species and
 *    reactions with low reaction flow rates, reducing the computational
 *    cost of the simulation.
 *
 *  - **DefinedEngineView:** An engine view that restricts the reaction
 *    network to a predefined set of species and reactions. This can be
 *    useful for simulating specific reaction pathways or for comparing
 *    results with other codes.
 *
 *  - **MultiscalePartitioningEngineView:** An engine view that partitions the
 *    reaction network into multiple groups based on timescales. This can be
 *    useful for simulating stiff reaction networks, where some reactions
 *    occur much faster than others.
 *
 *  - **NetworkPrimingEngineView:** An engine view that primes the reaction
 *    network with a specific species or set of species. This can be useful
 *    for igniting a reaction network or for studying the effects of specific
 *    species on the network.
 *
 * @section UsageExamples Usage Examples
 *
 * @subsection GraphEngineExample GraphEngine Example
 *
 * The following code shows how to create a GraphEngine from a composition:
 *
 * @code
 * #include "gridfire/engine/engine_graph.h"
 * #include "fourdst/composition/composition.h"
 *
 * // Create a composition
 * fourdst::composition::Composition composition;
 *
 * // Create a GraphEngine
 * gridfire::GraphEngine engine(composition);
 * @endcode
 *
 * @subsection AdaptiveEngineViewExample AdaptiveEngineView Example
 *
 * The following code shows how to create an AdaptiveEngineView from a
 * GraphEngine:
 *
 * @code
 * #include "gridfire/engine/views/engine_adaptive.h"
 * #include "gridfire/engine/engine_graph.h"
 * #include "fourdst/composition/composition.h"
 *
 * // Create a composition
 * fourdst::composition::Composition composition;
 *
 * // Create a GraphEngine
 * gridfire::GraphEngine baseEngine(composition);
 *
 * // Create an AdaptiveEngineView
 * gridfire::AdaptiveEngineView engine(baseEngine);
 * @endcode
 *
 * @subsection DefinedEngineViewExample DefinedEngineView Example
 *
 * The following code shows how to create a DefinedEngineView from a
 * GraphEngine:
 *
 * @code
 * #include "gridfire/engine/views/engine_defined.h"
 * #include "gridfire/engine/engine_graph.h"
 * #include "fourdst/composition/composition.h"
 *
 * // Create a composition
 * fourdst::composition::Composition composition;
 *
 * // Create a GraphEngine
 * gridfire::GraphEngine baseEngine(composition);
 *
 * // Create a DefinedEngineView
 * std::vector<std::string> peNames = {"p(p,e+)d", "he4(a,g)be8"};
 * gridfire::DefinedEngineView engine(peNames, baseEngine);
 * @endcode
 *
 * @subsection MultiscalePartitioningEngineViewExample MultiscalePartitioningEngineView Example
 *
 * The following code shows how to create a MultiscalePartitioningEngineView from a
 * GraphEngine:
 *
 * @code
 * #include "gridfire/engine/views/engine_multiscale.h"
 * #include "gridfire/engine/engine_graph.h"
 * #include "fourdst/composition/composition.h"
 *
 * // Create a composition
 * fourdst::composition::Composition composition;
 *
 * // Create a GraphEngine
 * gridfire::GraphEngine baseEngine(composition);
 *
 * // Create a MultiscalePartitioningEngineView
 * gridfire::MultiscalePartitioningEngineView engine(baseEngine);
 * @endcode
 *
 * @subsection NetworkPrimingEngineViewExample NetworkPrimingEngineView Example
 *
 * The following code shows how to create a NetworkPrimingEngineView from a
 * GraphEngine:
 *
 * @code
 * #include "gridfire/engine/views/engine_priming.h"
 * #include "gridfire/engine/engine_graph.h"
 * #include "fourdst/composition/composition.h"
 *
 * // Create a composition
 * fourdst::composition::Composition composition;
 *
 * // Create a GraphEngine
 * gridfire::GraphEngine baseEngine(composition);
 *
 * // Create a NetworkPrimingEngineView
 * std::string primingSymbol = "p";
 * gridfire::NetworkPrimingEngineView engine(primingSymbol, baseEngine);
 * @endcode
 */
#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/engine_graph.h"

#include "gridfire/engine/views/engine_views.h"
#include "gridfire/engine/procedures/engine_procedures.h"
#include "gridfire/engine/types/engine_types.h"