#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc if needed directly

#include <iostream>

#include "bindings.h"

#include "gridfire/engine/engine.h"
#include "trampoline/py_engine.h"



namespace py = pybind11;

namespace {
    template <typename T>
    concept IsDynamicEngine = std::is_base_of_v<gridfire::DynamicEngine, T>;

    template <IsDynamicEngine T, IsDynamicEngine BaseT>
    void registerDynamicEngineDefs(py::class_<T, BaseT> pyClass) {
        pyClass.def("calculateRHSAndEnergy", &T::calculateRHSAndEnergy,
            py::arg("Y"),
            py::arg("T9"),
            py::arg("rho"),
            "Calculate the right-hand side (dY/dt) and energy generation rate."
        )
        .def("generateJacobianMatrix", py::overload_cast<const std::vector<double>&, double, double>(&T::generateJacobianMatrix, py::const_),
            py::arg("Y_dynamic"),
            py::arg("T9"),
            py::arg("rho"),
            "Generate the Jacobian matrix for the current state."
        )
        .def("generateStoichiometryMatrix", &T::generateStoichiometryMatrix)
        .def("calculateMolarReactionFlow",
            static_cast<double (T::*)(const gridfire::reaction::Reaction&, const std::vector<double>&, const double, const double) const>(&T::calculateMolarReactionFlow),
            py::arg("reaction"),
            py::arg("Y"),
            py::arg("T9"),
            py::arg("rho"),
            "Calculate the molar reaction flow for a given reaction."
        )
        .def("getNetworkSpecies", &T::getNetworkSpecies,
            "Get the list of species in the network."
        )
        .def("getNetworkReactions", &T::getNetworkReactions,
            "Get the set of logical reactions in the network."
        )
        .def ("setNetworkReactions", &T::setNetworkReactions,
            py::arg("reactions"),
            "Set the network reactions to a new set of reactions."
        )
        .def("getJacobianMatrixEntry", &T::getJacobianMatrixEntry,
            py::arg("i"),
            py::arg("j"),
            "Get an entry from the previously generated Jacobian matrix."
        )
        .def("getStoichiometryMatrixEntry", &T::getStoichiometryMatrixEntry,
            py::arg("speciesIndex"),
            py::arg("reactionIndex"),
            "Get an entry from the stoichiometry matrix."
        )
        .def("getSpeciesTimescales", &T::getSpeciesTimescales,
            py::arg("Y"),
            py::arg("T9"),
            py::arg("rho"),
            "Get the timescales for each species in the network."
        )
        .def("getSpeciesDestructionTimescales", &T::getSpeciesDestructionTimescales,
            py::arg("Y"),
            py::arg("T9"),
            py::arg("rho"),
            "Get the destruction timescales for each species in the network."
        )
        .def("update", &T::update,
            py::arg("netIn"),
            "Update the engine state based on the provided NetIn object."
        )
        .def("setScreeningModel", &T::setScreeningModel,
            py::arg("screeningModel"),
            "Set the screening model for the engine."
        )
        .def("getScreeningModel", &T::getScreeningModel,
            "Get the current screening model of the engine."
        )
        .def("getSpeciesIndex", &T::getSpeciesIndex,
            py::arg("species"),
            "Get the index of a species in the network."
        )
        .def("mapNetInToMolarAbundanceVector", &T::mapNetInToMolarAbundanceVector,
            py::arg("netIn"),
            "Map a NetIn object to a vector of molar abundances."
        )
        .def("primeEngine", &T::primeEngine,
            py::arg("netIn"),
            "Prime the engine with a NetIn object to prepare for calculations."
        )
        .def("getDepth", &T::getDepth,
            "Get the current build depth of the engine."
        )
        .def("rebuild", &T::rebuild,
            py::arg("composition"),
            py::arg("depth") = gridfire::NetworkBuildDepth::Full,
            "Rebuild the engine with a new composition and build depth."
        )
        .def("isStale", &T::isStale,
            py::arg("netIn"),
            "Check if the engine is stale based on the provided NetIn object."
        );

    }
}

void register_engine_bindings(py::module &m) {
    register_base_engine_bindings(m);
    register_engine_view_bindings(m);

    m.def("build_reaclib_nuclear_network", &gridfire::build_reaclib_nuclear_network,
        py::arg("composition"),
        py::arg("maxLayers") = gridfire::NetworkBuildDepth::Full,
        py::arg("reverse") = false,
        "Build a nuclear network from a composition using ReacLib data."
    );

    py::enum_<gridfire::PrimingReportStatus>(m, "PrimingReportStatus")
        .value("FULL_SUCCESS", gridfire::PrimingReportStatus::FULL_SUCCESS, "Priming was full successful.")
        .value("NO_SPECIES_TO_PRIME", gridfire::PrimingReportStatus::NO_SPECIES_TO_PRIME, "No species to prime.")
        .value("MAX_ITERATIONS_REACHED", gridfire::PrimingReportStatus::MAX_ITERATIONS_REACHED, "Maximum iterations reached during priming.")
        .value("FAILED_TO_FINALIZE_COMPOSITION", gridfire::PrimingReportStatus::FAILED_TO_FINALIZE_COMPOSITION, "Failed to finalize the composition after priming.")
        .value("FAILED_TO_FIND_CREATION_CHANNEL", gridfire::PrimingReportStatus::FAILED_TO_FIND_CREATION_CHANNEL, "Failed to find a creation channel for the priming species.")
        .value("FAILED_TO_FIND_PRIMING_REACTIONS", gridfire::PrimingReportStatus::FAILED_TO_FIND_PRIMING_REACTIONS, "Failed to find priming reactions for the species.")
        .value("BASE_NETWORK_TOO_SHALLOW", gridfire::PrimingReportStatus::BASE_NETWORK_TOO_SHALLOW, "The base network is too shallow for priming.")
        .export_values()
        .def("__repr__", [](const gridfire::PrimingReportStatus& status) {
                std::stringstream ss;
                ss << gridfire::PrimingReportStatusStrings.at(status) << "\n";
                return ss.str();
            },
        "String representation of the PrimingReport."
        );

    py::class_<gridfire::PrimingReport>(m, "PrimingReport")
        .def_readonly("success", &gridfire::PrimingReport::success, "Indicates if the priming was successful.")
        .def_readonly("massFractionChanges", &gridfire::PrimingReport::massFractionChanges, "Map of species to their mass fraction changes after priming.")
        .def_readonly("primedComposition", &gridfire::PrimingReport::primedComposition, "The composition after priming.")
        .def_readonly("status", &gridfire::PrimingReport::status, "Status message from the priming process.")
        .def("__repr__", [](const gridfire::PrimingReport& report) {
                std::stringstream ss;
                ss << report;
                return ss.str();
            }
        );
}

void register_base_engine_bindings(pybind11::module &m) {

    py::class_<gridfire::StepDerivatives<double>>(m, "StepDerivatives")
        .def_readonly("dYdt", &gridfire::StepDerivatives<double>::dydt, "The right-hand side (dY/dt) of the ODE system.")
        .def_readonly("energy", &gridfire::StepDerivatives<double>::nuclearEnergyGenerationRate, "The energy generation rate.");

    py::class_<gridfire::SparsityPattern>(m, "SparsityPattern");

    abs_stype_register_engine_bindings(m);
    abs_stype_register_dynamic_engine_bindings(m);
    con_stype_register_graph_engine_bindings(m);
}

void abs_stype_register_engine_bindings(pybind11::module &m) {
    py::class_<gridfire::Engine, PyEngine>(m, "Engine");
}

void abs_stype_register_dynamic_engine_bindings(pybind11::module &m) {
    const auto a = py::class_<gridfire::DynamicEngine, PyDynamicEngine>(m, "DynamicEngine");
}

void con_stype_register_graph_engine_bindings(pybind11::module &m) {
    py::enum_<gridfire::NetworkBuildDepth>(m, "NetworkBuildDepth")
        .value("Full", gridfire::NetworkBuildDepth::Full, "Full network build depth")
        .value("Shallow", gridfire::NetworkBuildDepth::Shallow, "Shallow network build depth")
        .value("SecondOrder", gridfire::NetworkBuildDepth::SecondOrder, "Second order network build depth")
        .value("ThirdOrder", gridfire::NetworkBuildDepth::ThirdOrder, "Third order network build depth")
        .value("FourthOrder", gridfire::NetworkBuildDepth::FourthOrder, "Fourth order network build depth")
        .value("FifthOrder", gridfire::NetworkBuildDepth::FifthOrder, "Fifth order network build depth")
        .export_values();

    py::class_<gridfire::BuildDepthType>(m, "BuildDepthType");

    auto py_dynamic_engine_bindings = py::class_<gridfire::GraphEngine, gridfire::DynamicEngine>(m, "GraphEngine");

    // Register the Graph Engine Specific Bindings
    py_dynamic_engine_bindings.def(py::init<const fourdst::composition::Composition &, const gridfire::BuildDepthType>(),
        py::arg("composition"),
        py::arg("depth") = gridfire::NetworkBuildDepth::Full,
        "Initialize GraphEngine with a composition and build depth."
    );
    py_dynamic_engine_bindings.def(py::init<const fourdst::composition::Composition &, const gridfire::partition::PartitionFunction &, const gridfire::BuildDepthType>(),
        py::arg("composition"),
        py::arg("partitionFunction"),
        py::arg("depth") = gridfire::NetworkBuildDepth::Full,
        "Initialize GraphEngine with a composition, partition function and build depth."
    );
    py_dynamic_engine_bindings.def(py::init<const gridfire::reaction::LogicalReactionSet &>(),
        py::arg("reactions"),
        "Initialize GraphEngine with a set of reactions."
    );
    py_dynamic_engine_bindings.def("generateJacobianMatrix", py::overload_cast<const std::vector<double>&, double, double, const gridfire::SparsityPattern&>(&gridfire::GraphEngine::generateJacobianMatrix, py::const_),
        py::arg("Y_dynamic"),
        py::arg("T9"),
        py::arg("rho"),
        py::arg("sparsityPattern"),
        "Generate the Jacobian matrix for the current state with a specified sparsity pattern."
    );
    py_dynamic_engine_bindings.def_static("getNetReactionStoichiometry", &gridfire::GraphEngine::getNetReactionStoichiometry,
        py::arg("reaction"),
        "Get the net stoichiometry for a given reaction."
    );
    py_dynamic_engine_bindings.def("involvesSpecies", &gridfire::GraphEngine::involvesSpecies,
        py::arg("species"),
        "Check if a given species is involved in the network."
    );
    py_dynamic_engine_bindings.def("exportToDot", &gridfire::GraphEngine::exportToDot,
        py::arg("filename"),
        "Export the network to a DOT file for visualization."
    );
    py_dynamic_engine_bindings.def("exportToCSV", &gridfire::GraphEngine::exportToCSV,
        py::arg("filename"),
        "Export the network to a CSV file for analysis."
    );
    py_dynamic_engine_bindings.def("setPrecomputation", &gridfire::GraphEngine::setPrecomputation,
        py::arg("precompute"),
        "Enable or disable precomputation for the engine."
    );
    py_dynamic_engine_bindings.def("isPrecomputationEnabled", &gridfire::GraphEngine::isPrecomputationEnabled,
        "Check if precomputation is enabled for the engine."
    );
    py_dynamic_engine_bindings.def("getPartitionFunction", &gridfire::GraphEngine::getPartitionFunction,
        "Get the partition function used by the engine."
    );
    py_dynamic_engine_bindings.def("calculateReverseRate", &gridfire::GraphEngine::calculateReverseRate,
        py::arg("reaction"),
        py::arg("T9"),
        "Calculate the reverse rate for a given reaction at a specific temperature."
    );
    py_dynamic_engine_bindings.def("calculateReverseRateTwoBody", &gridfire::GraphEngine::calculateReverseRateTwoBody,
        py::arg("reaction"),
        py::arg("T9"),
        py::arg("forwardRate"),
        py::arg("expFactor"),
        "Calculate the reverse rate for a two-body reaction at a specific temperature."
    );
    py_dynamic_engine_bindings.def("calculateReverseRateTwoBodyDerivative", &gridfire::GraphEngine::calculateReverseRateTwoBodyDerivative,
        py::arg("reaction"),
        py::arg("T9"),
        py::arg("reverseRate"),
        "Calculate the derivative of the reverse rate for a two-body reaction at a specific temperature."
    );
    py_dynamic_engine_bindings.def("isUsingReverseReactions", &gridfire::GraphEngine::isUsingReverseReactions,
        "Check if the engine is using reverse reactions."
    );
    py_dynamic_engine_bindings.def("setUseReverseReactions", &gridfire::GraphEngine::setUseReverseReactions,
        py::arg("useReverse"),
        "Enable or disable the use of reverse reactions in the engine."
    );


    // Register the general dynamic engine bindings
    registerDynamicEngineDefs<gridfire::GraphEngine, gridfire::DynamicEngine>(py_dynamic_engine_bindings);
}

void register_engine_view_bindings(pybind11::module &m) {
    auto py_defined_engine_view_bindings = py::class_<gridfire::DefinedEngineView, gridfire::DynamicEngine>(m, "DefinedEngineView");

    py_defined_engine_view_bindings.def(py::init<std::vector<std::string>, gridfire::DynamicEngine&>(),
        py::arg("peNames"),
        py::arg("baseEngine"),
        "Construct a defined engine view with a list of tracked reactions and a base engine.");
    py_defined_engine_view_bindings.def("getBaseEngine", &gridfire::DefinedEngineView::getBaseEngine,
        "Get the base engine associated with this defined engine view.");

    registerDynamicEngineDefs<gridfire::DefinedEngineView, gridfire::DynamicEngine>(py_defined_engine_view_bindings);

    auto py_file_defined_engine_view_bindings = py::class_<gridfire::FileDefinedEngineView, gridfire::DefinedEngineView>(m, "FileDefinedEngineView");
    py_file_defined_engine_view_bindings.def(py::init<gridfire::DynamicEngine&, const std::string&, const gridfire::io::NetworkFileParser&>(),
        py::arg("baseEngine"),
        py::arg("fileName"),
        py::arg("parser"),
        "Construct a defined engine view from a file and a base engine."
    );
    py_file_defined_engine_view_bindings.def("getNetworkFile", &gridfire::FileDefinedEngineView::getNetworkFile,
        "Get the network file associated with this defined engine view."
    );
    py_file_defined_engine_view_bindings.def("getParser", &gridfire::FileDefinedEngineView::getParser,
        "Get the parser used for this defined engine view."
    );
    py_file_defined_engine_view_bindings.def("getBaseEngine", &gridfire::FileDefinedEngineView::getBaseEngine,
        "Get the base engine associated with this file defined engine view.");

    registerDynamicEngineDefs<gridfire::FileDefinedEngineView, gridfire::DefinedEngineView>(py_file_defined_engine_view_bindings);

    auto py_priming_engine_view_bindings = py::class_<gridfire::NetworkPrimingEngineView, gridfire::DefinedEngineView>(m, "NetworkPrimingEngineView");
    py_priming_engine_view_bindings.def(py::init<const std::string&, gridfire::DynamicEngine&>(),
        py::arg("primingSymbol"),
        py::arg("baseEngine"),
        "Construct a priming engine view with a priming symbol and a base engine.");
    py_priming_engine_view_bindings.def(py::init<const fourdst::atomic::Species&, gridfire::DynamicEngine&>(),
        py::arg("primingSpecies"),
        py::arg("baseEngine"),
        "Construct a priming engine view with a priming species and a base engine.");
    py_priming_engine_view_bindings.def("getBaseEngine", &gridfire::NetworkPrimingEngineView::getBaseEngine,
        "Get the base engine associated with this priming engine view.");

    registerDynamicEngineDefs<gridfire::NetworkPrimingEngineView, gridfire::DefinedEngineView>(py_priming_engine_view_bindings);

    auto py_adaptive_engine_view_bindings = py::class_<gridfire::AdaptiveEngineView, gridfire::DynamicEngine>(m, "AdaptiveEngineView");
    py_adaptive_engine_view_bindings.def(py::init<gridfire::DynamicEngine&>(),
        py::arg("baseEngine"),
        "Construct an adaptive engine view with a base engine.");
    py_adaptive_engine_view_bindings.def("getBaseEngine", &gridfire::AdaptiveEngineView::getBaseEngine,
        "Get the base engine associated with this adaptive engine view.");

    registerDynamicEngineDefs<gridfire::AdaptiveEngineView, gridfire::DynamicEngine>(py_adaptive_engine_view_bindings);

    auto py_qse_cache_config = py::class_<gridfire::QSECacheConfig>(m, "QSECacheConfig");
    auto py_qse_cache_key = py::class_<gridfire::QSECacheKey>(m, "QSECacheKey");

    py_qse_cache_key.def(py::init<double, double, const std::vector<double>&>(),
        py::arg("T9"),
        py::arg("rho"),
        py::arg("Y")
    );

    py_qse_cache_key.def("hash", &gridfire::QSECacheKey::hash,
        "Get the pre-computed hash value of the key");

    py_qse_cache_key.def_static("bin", &gridfire::QSECacheKey::bin,
        py::arg("value"),
        py::arg("tol"),
        "bin a value based on a tolerance");
    py_qse_cache_key.def("__eq__", &gridfire::QSECacheKey::operator==,
        py::arg("other"),
        "Check if two QSECacheKeys are equal");

    auto py_multiscale_engine_view_bindings = py::class_<gridfire::MultiscalePartitioningEngineView, gridfire::DynamicEngine>(m, "MultiscalePartitioningEngineView");
    py_multiscale_engine_view_bindings.def(py::init<gridfire::GraphEngine&>(),
        py::arg("baseEngine"),
        "Construct a multiscale partitioning engine view with a base engine.");
    py_multiscale_engine_view_bindings.def("getBaseEngine", &gridfire::MultiscalePartitioningEngineView::getBaseEngine,
        "Get the base engine associated with this multiscale partitioning engine view.");
    py_multiscale_engine_view_bindings.def("analyzeTimescalePoolConnectivity", &gridfire::MultiscalePartitioningEngineView::analyzeTimescalePoolConnectivity,
        py::arg("timescale_pools"),
        py::arg("Y"),
        py::arg("T9"),
        py::arg("rho"),
        "Analyze the connectivity of timescale pools in the network.");
    py_multiscale_engine_view_bindings.def("partitionNetwork", py::overload_cast<const std::vector<double>&, double, double>(&gridfire::MultiscalePartitioningEngineView::partitionNetwork),
        py::arg("Y"),
        py::arg("T9"),
        py::arg("rho"),
        "Partition the network based on species timescales and connectivity.");
    py_multiscale_engine_view_bindings.def("partitionNetwork", py::overload_cast<const gridfire::NetIn&>(&gridfire::MultiscalePartitioningEngineView::partitionNetwork),
        py::arg("netIn"),
        "Partition the network based on a NetIn object.");
    py_multiscale_engine_view_bindings.def("exportToDot", &gridfire::MultiscalePartitioningEngineView::exportToDot,
        py::arg("filename"),
        py::arg("Y"),
        py::arg("T9"),
        py::arg("rho"),
        "Export the network to a DOT file for visualization.");
    py_multiscale_engine_view_bindings.def("getFastSpecies", &gridfire::MultiscalePartitioningEngineView::getFastSpecies,
        "Get the list of fast species in the network.");
    py_multiscale_engine_view_bindings.def("getDynamicSpecies", &gridfire::MultiscalePartitioningEngineView::getDynamicSpecies,
        "Get the list of dynamic species in the network.");
    py_multiscale_engine_view_bindings.def("equilibrateNetwork", py::overload_cast<const std::vector<double>&, double, double>(&gridfire::MultiscalePartitioningEngineView::equilibrateNetwork),
        py::arg("Y"),
        py::arg("T9"),
        py::arg("rho"),
        "Equilibrate the network based on species abundances and conditions.");
    py_multiscale_engine_view_bindings.def("equilibrateNetwork", py::overload_cast<const gridfire::NetIn&>(&gridfire::MultiscalePartitioningEngineView::equilibrateNetwork),
        py::arg("netIn"),
        "Equilibrate the network based on a NetIn object.");

    registerDynamicEngineDefs<gridfire::MultiscalePartitioningEngineView, gridfire::DynamicEngine>(py_multiscale_engine_view_bindings);
}








