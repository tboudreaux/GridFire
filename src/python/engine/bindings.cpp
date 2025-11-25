#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly

#include "bindings.h"

#include "gridfire/engine/engine.h"
#include "gridfire/exceptions/exceptions.h"
#include "pybind11/numpy.h"
#include "trampoline/py_engine.h"



namespace py = pybind11;

namespace {
    template <typename T>
    concept IsDynamicEngine = std::is_base_of_v<gridfire::engine::DynamicEngine, T>;

    template <IsDynamicEngine T, IsDynamicEngine BaseT>
    void registerDynamicEngineDefs(py::class_<T, BaseT> pyClass) {
        pyClass.def(
            "calculateRHSAndEnergy",
            [](
                const gridfire::engine::DynamicEngine& self,
                const fourdst::composition::Composition& comp,
                const double T9,
                const double rho
            ) {
                auto result = self.calculateRHSAndEnergy(comp, T9, rho);
                if (!result.has_value()) {
                    throw gridfire::exceptions::EngineError(std::format("calculateRHSAndEnergy returned a potentially recoverable error {}", gridfire::engine::EngineStatus_to_string(result.error())));
                }
                return result.value();
            },
            py::arg("comp"),
            py::arg("T9"),
            py::arg("rho"),
            "Calculate the right-hand side (dY/dt) and energy generation rate."
        )
        .def("calculateEpsDerivatives",
            &gridfire::engine::DynamicEngine::calculateEpsDerivatives,
            py::arg("comp"),
            py::arg("T9"),
            py::arg("rho"),
            "Calculate deps/dT and deps/drho"
        )
        .def("generateJacobianMatrix",
            [](const gridfire::engine::DynamicEngine& self,
               const fourdst::composition::Composition& comp,
               const double T9,
               const double rho) -> gridfire::engine::NetworkJacobian {
                return self.generateJacobianMatrix(comp, T9, rho);
            },
            py::arg("comp"),
            py::arg("T9"),
            py::arg("rho"),
            "Generate the Jacobian matrix for the current state."
        )
        .def("generateJacobianMatrix",
            [](const gridfire::engine::DynamicEngine& self,
               const fourdst::composition::Composition& comp,
               const double T9,
               const double rho,
               const std::vector<fourdst::atomic::Species>& activeSpecies) -> gridfire::engine::NetworkJacobian {
                return self.generateJacobianMatrix(comp, T9, rho, activeSpecies);
            },
            py::arg("comp"),
            py::arg("T9"),
            py::arg("rho"),
            py::arg("activeSpecies"),
            "Generate the jacobian matrix only for the subset of the matrix representing the active species."
        )
        .def("generateJacobianMatrix",
            [](const gridfire::engine::DynamicEngine& self,
               const fourdst::composition::Composition& comp,
               const double T9,
               const double rho,
               const gridfire::engine::SparsityPattern& sparsityPattern) -> gridfire::engine::NetworkJacobian {
                return self.generateJacobianMatrix(comp, T9, rho, sparsityPattern);
            },
            py::arg("comp"),
            py::arg("T9"),
            py::arg("rho"),
            py::arg("sparsityPattern"),
            "Generate the jacobian matrix for the given sparsity pattern"
        )
        .def("generateStoichiometryMatrix",
            &T::generateStoichiometryMatrix
        )
        .def("calculateMolarReactionFlow",
            [](
                const gridfire::engine::DynamicEngine& self,
                const gridfire::reaction::Reaction& reaction,
                const fourdst::composition::Composition& comp,
                const double T9,
                const double rho
            ) -> double {
                return self.calculateMolarReactionFlow(reaction, comp, T9, rho);
            },
            py::arg("reaction"),
            py::arg("comp"),
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
        .def("getStoichiometryMatrixEntry", &T::getStoichiometryMatrixEntry,
            py::arg("species"),
            py::arg("reaction"),
            "Get an entry from the stoichiometry matrix."
        )
        .def("getSpeciesTimescales",
            [](
                const gridfire::engine::DynamicEngine& self,
                const fourdst::composition::Composition& comp,
                const double T9,
                const double rho
            ) -> std::unordered_map<fourdst::atomic::Species, double> {
                const auto result = self.getSpeciesTimescales(comp, T9, rho);
                if (!result.has_value()) {
                    throw gridfire::exceptions::EngineError(std::format("getSpeciesTimescales has returned a potentially recoverable error {}", gridfire::engine::EngineStatus_to_string(result.error())));
                }
                return result.value();
            },
            py::arg("comp"),
            py::arg("T9"),
            py::arg("rho"),
            "Get the timescales for each species in the network."
        )
        .def("getSpeciesDestructionTimescales",
            [](
                const gridfire::engine::DynamicEngine& self,
                const fourdst::composition::Composition& comp,
                const double T9,
                const double rho
            ) -> std::unordered_map<fourdst::atomic::Species, double> {
                const auto result = self.getSpeciesDestructionTimescales(comp, T9, rho);
                if (!result.has_value()) {
                    throw gridfire::exceptions::EngineError(std::format("getSpeciesDestructionTimescales has returned a potentially recoverable error {}", gridfire::engine::EngineStatus_to_string(result.error())));
                }
                return result.value();
            },
            py::arg("comp"),
            py::arg("T9"),
            py::arg("rho"),
            "Get the destruction timescales for each species in the network."
        )
        .def("update",
            &T::update,
            py::arg("netIn"),
            "Update the engine state based on the provided NetIn object."
        )
        .def("setScreeningModel",
            &T::setScreeningModel,
            py::arg("screeningModel"),
            "Set the screening model for the engine."
        )
        .def("getScreeningModel",
            &T::getScreeningModel,
            "Get the current screening model of the engine."
        )
        .def("getSpeciesIndex",
            &T::getSpeciesIndex,
            py::arg("species"),
            "Get the index of a species in the network."
        )
        .def("mapNetInToMolarAbundanceVector",
            &T::mapNetInToMolarAbundanceVector,
            py::arg("netIn"),
            "Map a NetIn object to a vector of molar abundances."
        )
        .def("primeEngine",
            &T::primeEngine,
            py::arg("netIn"),
            "Prime the engine with a NetIn object to prepare for calculations."
        )
        .def("getDepth",
            &T::getDepth,
            "Get the current build depth of the engine."
        )
        .def("rebuild",
            &T::rebuild,
            py::arg("composition"),
            py::arg("depth") = gridfire::engine::NetworkBuildDepth::Full,
            "Rebuild the engine with a new composition and build depth."
        )
        .def("isStale",
            &T::isStale,
            py::arg("netIn"),
            "Check if the engine is stale based on the provided NetIn object."
        )
        .def("collectComposition",
            &T::collectComposition,
            py::arg("composition"),
            py::arg("T9"),
            py::arg("rho"),
            "Recursively collect composition from current engine and any sub engines if they exist."
        )
        .def("getSpeciesStatus",
            &T::getSpeciesStatus,
            py::arg("species"),
            "Get the status of a species in the network."
        );

    }
}

void register_engine_bindings(py::module &m) {
    register_engine_type_bindings(m);
    register_engine_procedural_bindings(m);
    register_base_engine_bindings(m);
    register_engine_view_bindings(m);
    register_engine_diagnostic_bindings(m);
}

void register_base_engine_bindings(const pybind11::module &m) {

    py::class_<gridfire::engine::StepDerivatives<double>>(m, "StepDerivatives")
        .def_readonly("dYdt", &gridfire::engine::StepDerivatives<double>::dydt, "The right-hand side (dY/dt) of the ODE system.")
        .def_readonly("energy", &gridfire::engine::StepDerivatives<double>::nuclearEnergyGenerationRate, "The energy generation rate.");

    py::class_<gridfire::engine::SparsityPattern> py_sparsity_pattern(m, "SparsityPattern");

    abs_stype_register_engine_bindings(m);
    abs_stype_register_dynamic_engine_bindings(m);
    con_stype_register_graph_engine_bindings(m);
}

void abs_stype_register_engine_bindings(const pybind11::module &m) {
    py::class_<gridfire::engine::Engine, PyEngine>(m, "Engine");
}

void abs_stype_register_dynamic_engine_bindings(const pybind11::module &m) {
    const auto a = py::class_<gridfire::engine::DynamicEngine, PyDynamicEngine>(m, "DynamicEngine");
}

void register_engine_procedural_bindings(pybind11::module &m) {
    register_engine_construction_bindings(m);
    register_engine_priming_bindings(m);
}

void register_engine_diagnostic_bindings(pybind11::module &m) {
    auto diagnostics = m.def_submodule("diagnostics", "A submodule for engine diagnostics");
    diagnostics.def("report_limiting_species",
        &gridfire::engine::diagnostics::report_limiting_species,
        py::arg("engine"),
        py::arg("Y_full"),
        py::arg("E_full"),
        py::arg("relTol"),
        py::arg("absTol"),
        py::arg("top_n"),
        py::arg("json")
    );

    diagnostics.def("inspect_species_balance",
        &gridfire::engine::diagnostics::inspect_species_balance,
        py::arg("engine"),
        py::arg("species_name"),
        py::arg("comp"),
        py::arg("T9"),
        py::arg("rho"),
        py::arg("json")
    );

    diagnostics.def("inspect_jacobian_stiffness",
        &gridfire::engine::diagnostics::inspect_jacobian_stiffness,
        py::arg("engine"),
        py::arg("comp"),
        py::arg("T9"),
        py::arg("rho"),
        py::arg("json")
    );
}

void register_engine_construction_bindings(pybind11::module &m) {
    py::enum_<gridfire::engine::NetworkConstructionFlags>(m, "NetworkConstructionFlags")
        .value("NONE", gridfire::engine::NetworkConstructionFlags::NONE, "No special construction flags.")
        .value("REACLIB_STRONG", gridfire::engine::NetworkConstructionFlags::REACLIB_STRONG, "Include strong reactions from reaclib.")
        .value("WRL_BETA_MINUS", gridfire::engine::NetworkConstructionFlags::WRL_BETA_MINUS, "Include beta-minus decay reactions from weak rate library.")
        .value("WRL_BETA_PLUS", gridfire::engine::NetworkConstructionFlags::WRL_BETA_PLUS, "Include beta-plus decay reactions from weak rate library.")
        .value("WRL_ELECTRON_CAPTURE", gridfire::engine::NetworkConstructionFlags::WRL_ELECTRON_CAPTURE, "Include electron capture reactions from weak rate library.")
        .value("WRL_POSITRON_CAPTURE", gridfire::engine::NetworkConstructionFlags::WRL_POSITRON_CAPTURE, "Include positron capture reactions from weak rate library.")
        .value("REACLIB_WEAK", gridfire::engine::NetworkConstructionFlags::REACLIB_WEAK, "Include weak reactions from reaclib.")
        .value("WRL_WEAK", gridfire::engine::NetworkConstructionFlags::WRL_WEAK, "Include all weak reactions from weak rate library.")
        .value("REACLIB", gridfire::engine::NetworkConstructionFlags::REACLIB, "Include all reactions from reaclib.")
        .value("DEFAULT", gridfire::engine::NetworkConstructionFlags::DEFAULT, "Default construction flags (Reaclib strong and weak).")
        .export_values()
        .def("__repr__", [](const gridfire::engine::NetworkConstructionFlags& flags) {
            return gridfire::engine::NetworkConstructionFlagsToString(flags);
        });

    m.def("build_nuclear_network", &gridfire::engine::build_nuclear_network,
        py::arg("composition"),
        py::arg("weakInterpolator"),
        py::arg("maxLayers") = gridfire::engine::NetworkBuildDepth::Full,
        py::arg("ReactionTypes") = gridfire::engine::NetworkConstructionFlags::DEFAULT,
        "Build a nuclear network from a composition using all archived reaction data."
    );
}

void register_engine_priming_bindings(pybind11::module &m) {
    m.def("primeNetwork",
        &gridfire::engine::primeNetwork,
        py::arg("netIn"),
        py::arg("engine"),
        py::arg("ignoredReactionTypes") = std::nullopt,
        "Prime a network with a short timescale ignition"
    );
}

void register_engine_type_bindings(pybind11::module &m) {
    register_engine_building_type_bindings(m);
    register_engine_reporting_type_bindings(m);
    register_engine_types_bindings(m);
    register_jacobian_type_bindings(m);

}

void register_engine_building_type_bindings(const pybind11::module &m) {
    py::enum_<gridfire::engine::NetworkBuildDepth>(m, "NetworkBuildDepth")
        .value("Full", gridfire::engine::NetworkBuildDepth::Full, "Full network build depth")
        .value("Shallow", gridfire::engine::NetworkBuildDepth::Shallow, "Shallow network build depth")
        .value("SecondOrder", gridfire::engine::NetworkBuildDepth::SecondOrder, "Second order network build depth")
        .value("ThirdOrder", gridfire::engine::NetworkBuildDepth::ThirdOrder, "Third order network build depth")
        .value("FourthOrder", gridfire::engine::NetworkBuildDepth::FourthOrder, "Fourth order network build depth")
        .value("FifthOrder", gridfire::engine::NetworkBuildDepth::FifthOrder, "Fifth order network build depth")
        .export_values();

    py::class_<gridfire::engine::BuildDepthType> py_build_depth_type(m, "BuildDepthType");
}

void register_engine_reporting_type_bindings(const pybind11::module &m) {
    py::enum_<gridfire::engine::PrimingReportStatus>(m, "PrimingReportStatus")
        .value("FULL_SUCCESS", gridfire::engine::PrimingReportStatus::SUCCESS, "Priming was full successful.")
        .value("NO_SPECIES_TO_PRIME", gridfire::engine::PrimingReportStatus::SOLVER_FAILURE, "Solver Failed to converge during priming.")
        .value("MAX_ITERATIONS_REACHED", gridfire::engine::PrimingReportStatus::ALREADY_PRIMED, "Engine has already been primed.")
        .export_values()
        .def("__repr__", [](const gridfire::engine::PrimingReportStatus& status) {
                std::stringstream ss;
                ss << gridfire::engine::PrimingReportStatusStrings.at(status) << "\n";
                return ss.str();
            },
        "String representation of the PrimingReport."
        );

    py::class_<gridfire::engine::PrimingReport>(m, "PrimingReport")
        .def_readonly("success", &gridfire::engine::PrimingReport::success, "Indicates if the priming was successful.")
        .def_readonly("primedComposition", &gridfire::engine::PrimingReport::primedComposition, "The composition after priming.")
        .def_readonly("status", &gridfire::engine::PrimingReport::status, "Status message from the priming process.")
        .def("__repr__", [](const gridfire::engine::PrimingReport& report) {
                std::stringstream ss;
                ss << report;
                return ss.str();
            }
        );

    py::enum_<gridfire::engine::SpeciesStatus>(m, "SpeciesStatus")
        .value("ACTIVE", gridfire::engine::SpeciesStatus::ACTIVE, "Species is active in the network.")
        .value("EQUILIBRIUM", gridfire::engine::SpeciesStatus::EQUILIBRIUM, "Species is in equilibrium.")
        .value("INACTIVE_FLOW", gridfire::engine::SpeciesStatus::INACTIVE_FLOW, "Species is inactive due to flow.")
        .value("NOT_PRESENT", gridfire::engine::SpeciesStatus::NOT_PRESENT, "Species is not present in the network.")
        .export_values()
        .def("__repr__", [](const gridfire::engine::SpeciesStatus& status) {
            return gridfire::engine::SpeciesStatus_to_string(status);
        });
}

void register_engine_types_bindings(const pybind11::module &m) {
    py::enum_<gridfire::engine::EngineTypes>(m, "EngineTypes")
        .value("GRAPH_ENGINE", gridfire::engine::EngineTypes::GRAPH_ENGINE, "The standard graph-based engine.")
        .value("ADAPTIVE_ENGINE_VIEW", gridfire::engine::EngineTypes::ADAPTIVE_ENGINE_VIEW, "An engine that adapts based on certain criteria.")
        .value("MULTISCALE_PARTITIONING_ENGINE_VIEW", gridfire::engine::EngineTypes::MULTISCALE_PARTITIONING_ENGINE_VIEW, "An engine that partitions the system at multiple scales.")
        .value("PRIMING_ENGINE_VIEW", gridfire::engine::EngineTypes::PRIMING_ENGINE_VIEW, "An engine that uses a priming strategy for simulations.")
        .value("DEFINED_ENGINE_VIEW", gridfire::engine::EngineTypes::DEFINED_ENGINE_VIEW, "An engine defined by user specifications.")
        .value("FILE_DEFINED_ENGINE_VIEW", gridfire::engine::EngineTypes::FILE_DEFINED_ENGINE_VIEW, "An engine defined through external files.")
        .export_values()
        .def("__repr__", [](const gridfire::engine::EngineTypes& type) {
                return std::string(gridfire::engine::engine_type_to_string(type));
            },
        "String representation of the EngineTypes."
        );
}

void register_jacobian_type_bindings(pybind11::module &m) {
    py::class_<gridfire::engine::NetworkJacobian>(m, "NetworkJacobian")
        .def("rank", &gridfire::engine::NetworkJacobian::rank, "Get the rank of the Jacobian matrix.")
        .def("nnz", &gridfire::engine::NetworkJacobian::nnz, "Get the number of non-zero entries in the Jacobian matrix.")
        .def("singular", &gridfire::engine::NetworkJacobian::singular, "Check if the Jacobian matrix is singular.")
        .def("infs", &gridfire::engine::NetworkJacobian::infs, "Get all infinite entries in the Jacobian matrix.")
        .def("nans", &gridfire::engine::NetworkJacobian::nans, "Get all NaN entries in the Jacobian matrix.")
        .def("data", &gridfire::engine::NetworkJacobian::data, "Get the underlying sparse matrix data.")
        .def("mapping", &gridfire::engine::NetworkJacobian::mapping, "Get the species-to-index mapping.")
        .def("shape", &gridfire::engine::NetworkJacobian::shape, "Get the shape of the Jacobian matrix as (rows, columns).")
        .def("to_csv", &gridfire::engine::NetworkJacobian::to_csv, py::arg("filename"), "Export the Jacobian matrix to a CSV file.")
        .def("__getitem__", [](const gridfire::engine::NetworkJacobian &self, const std::pair<fourdst::atomic::Species, fourdst::atomic::Species>& speciesPair) {
            return self(speciesPair.first, speciesPair.second);
        }, py::arg("key"), "Get an entry from the Jacobian matrix using species identifiers.")
        .def("__getitem__", [](const gridfire::engine::NetworkJacobian &self, const std::pair<size_t, size_t>& indexPair) {
            return self(indexPair.first, indexPair.second);
        }, py::arg("key"), "Get an entry from the Jacobian matrix using indices.")
        .def("__setitem__", [](gridfire::engine::NetworkJacobian &self, const std::pair<fourdst::atomic::Species, fourdst::atomic::Species>& speciesPair, double value) {
            self.set(speciesPair.first, speciesPair.second, value);
        }, py::arg("key"), py::arg("value"), "Set an entry in the Jacobian matrix using species identifiers.")
        .def("__setitem__", [](gridfire::engine::NetworkJacobian &self, const std::pair<size_t, size_t>& indexPair, double value) {
            self.set(indexPair.first, indexPair.second, value);
        }, py::arg("key"), py::arg("value"), "Set an entry in the Jacobian matrix using indices.")
        .def("to_numpy", [](const gridfire::engine::NetworkJacobian &self) {
            auto denseMatrix = Eigen::MatrixXd(self.data());
            return pybind11::array_t<double>(
                {std::get<0>(self.shape()), std::get<1>(self.shape())},
                {sizeof(double) * std::get<1>(self.shape()), sizeof(double)},
                denseMatrix.data()
            );
        }, "Convert the Jacobian matrix to a NumPy array.");

    m.def(
        "regularize_jacobian",
        [](const gridfire::engine::NetworkJacobian& jac, const fourdst::composition::Composition& comp) {
            return regularize_jacobian(jac, comp, std::nullopt);
        },
        py::arg("jacobian"),
        py::arg("composition"),
        "regularize_jacobian"
    );
}



void con_stype_register_graph_engine_bindings(const pybind11::module &m) {

    auto py_graph_engine_bindings = py::class_<gridfire::engine::GraphEngine, gridfire::engine::DynamicEngine>(m, "GraphEngine");

    // Register the Graph Engine Specific Bindings
    py_graph_engine_bindings.def(py::init<const fourdst::composition::Composition &, const gridfire::engine::BuildDepthType>(),
        py::arg("composition"),
        py::arg("depth") = gridfire::engine::NetworkBuildDepth::Full,
        "Initialize GraphEngine with a composition and build depth."
    );
    py_graph_engine_bindings.def(py::init<const fourdst::composition::Composition &,const gridfire::partition::PartitionFunction &, const gridfire::engine::BuildDepthType>(),
        py::arg("composition"),
        py::arg("partitionFunction"),
        py::arg("depth") = gridfire::engine::NetworkBuildDepth::Full,
        "Initialize GraphEngine with a composition, partition function and build depth."
    );
    py_graph_engine_bindings.def(py::init<const gridfire::reaction::ReactionSet &>(),
        py::arg("reactions"),
        "Initialize GraphEngine with a set of reactions."
    );
    py_graph_engine_bindings.def_static("getNetReactionStoichiometry",
        &gridfire::engine::GraphEngine::getNetReactionStoichiometry,
        py::arg("reaction"),
        "Get the net stoichiometry for a given reaction."
    );
    py_graph_engine_bindings.def("getSpeciesTimescales",
        [](const gridfire::engine::GraphEngine& self,
           const fourdst::composition::Composition& composition,
           const double T9,
           const double rho,
           const gridfire::reaction::ReactionSet& activeReactions) {
            return self.getSpeciesTimescales(composition, T9, rho, activeReactions);
        },
        py::arg("composition"),
        py::arg("T9"),
        py::arg("rho"),
        py::arg("activeReactions")
    );
    py_graph_engine_bindings.def("getSpeciesDestructionTimescales",
        [](const gridfire::engine::GraphEngine& self,
           const fourdst::composition::Composition& composition,
           const double T9,
           const double rho,
           const gridfire::reaction::ReactionSet& activeReactions) {
            return self.getSpeciesDestructionTimescales(composition, T9, rho, activeReactions);
        },
        py::arg("composition"),
        py::arg("T9"),
        py::arg("rho"),
        py::arg("activeReactions")
    );
    py_graph_engine_bindings.def("involvesSpecies",
        &gridfire::engine::GraphEngine::involvesSpecies,
        py::arg("species"),
        "Check if a given species is involved in the network."
    );
    py_graph_engine_bindings.def("exportToDot",
        &gridfire::engine::GraphEngine::exportToDot,
        py::arg("filename"),
        "Export the network to a DOT file for visualization."
    );
    py_graph_engine_bindings.def("exportToCSV",
        &gridfire::engine::GraphEngine::exportToCSV,
        py::arg("filename"),
        "Export the network to a CSV file for analysis."
    );
    py_graph_engine_bindings.def("setPrecomputation",
        &gridfire::engine::GraphEngine::setPrecomputation,
        py::arg("precompute"),
        "Enable or disable precomputation for the engine."
    );
    py_graph_engine_bindings.def("isPrecomputationEnabled",
        &gridfire::engine::GraphEngine::isPrecomputationEnabled,
        "Check if precomputation is enabled for the engine."
    );
    py_graph_engine_bindings.def("getPartitionFunction",
        &gridfire::engine::GraphEngine::getPartitionFunction,
        "Get the partition function used by the engine."
    );
    py_graph_engine_bindings.def("calculateReverseRate",
        &gridfire::engine::GraphEngine::calculateReverseRate,
        py::arg("reaction"),
        py::arg("T9"),
        py::arg("rho"),
        py::arg("composition"),
        "Calculate the reverse rate for a given reaction at a specific temperature, density, and composition."
    );
    py_graph_engine_bindings.def("calculateReverseRateTwoBody",
        &gridfire::engine::GraphEngine::calculateReverseRateTwoBody,
        py::arg("reaction"),
        py::arg("T9"),
        py::arg("forwardRate"),
        py::arg("expFactor"),
        "Calculate the reverse rate for a two-body reaction at a specific temperature."
    );
    py_graph_engine_bindings.def("calculateReverseRateTwoBodyDerivative",
        &gridfire::engine::GraphEngine::calculateReverseRateTwoBodyDerivative,
        py::arg("reaction"),
        py::arg("T9"),
        py::arg("rho"),
        py::arg("composition"),
        py::arg("reverseRate"),
        "Calculate the derivative of the reverse rate for a two-body reaction at a specific temperature."
    );
    py_graph_engine_bindings.def("isUsingReverseReactions",
        &gridfire::engine::GraphEngine::isUsingReverseReactions,
        "Check if the engine is using reverse reactions."
    );
    py_graph_engine_bindings.def("setUseReverseReactions",
        &gridfire::engine::GraphEngine::setUseReverseReactions,
        py::arg("useReverse"),
        "Enable or disable the use of reverse reactions in the engine."
    );

    // Register the general dynamic engine bindings
    registerDynamicEngineDefs<gridfire::engine::GraphEngine, gridfire::engine::DynamicEngine>(py_graph_engine_bindings);
}

void register_engine_view_bindings(const pybind11::module &m) {
    auto py_defined_engine_view_bindings = py::class_<gridfire::engine::DefinedEngineView, gridfire::engine::DynamicEngine>(m, "DefinedEngineView");

    py_defined_engine_view_bindings.def(py::init<std::vector<std::string>, gridfire::engine::GraphEngine&>(),
        py::arg("peNames"),
        py::arg("baseEngine"),
        "Construct a defined engine view with a list of tracked reactions and a base engine."
    );
    py_defined_engine_view_bindings.def("getBaseEngine", &gridfire::engine::DefinedEngineView::getBaseEngine,
        "Get the base engine associated with this defined engine view.");

    registerDynamicEngineDefs<gridfire::engine::DefinedEngineView, gridfire::engine::DynamicEngine>(py_defined_engine_view_bindings);

    auto py_file_defined_engine_view_bindings = py::class_<gridfire::engine::FileDefinedEngineView, gridfire::engine::DefinedEngineView>(m, "FileDefinedEngineView");
    py_file_defined_engine_view_bindings.def(
    py::init<gridfire::engine::GraphEngine&, const std::string&, const gridfire::io::NetworkFileParser&>(),
        py::arg("baseEngine"),
        py::arg("fileName"),
        py::arg("parser"),
        "Construct a defined engine view from a file and a base engine."
    );
    py_file_defined_engine_view_bindings.def("getNetworkFile", &gridfire::engine::FileDefinedEngineView::getNetworkFile,
        "Get the network file associated with this defined engine view."
    );
    py_file_defined_engine_view_bindings.def("getParser", &gridfire::engine::FileDefinedEngineView::getParser,
        "Get the parser used for this defined engine view."
    );
    py_file_defined_engine_view_bindings.def("getBaseEngine", &gridfire::engine::FileDefinedEngineView::getBaseEngine,
        "Get the base engine associated with this file defined engine view.");

    registerDynamicEngineDefs<gridfire::engine::FileDefinedEngineView, gridfire::engine::DefinedEngineView>(py_file_defined_engine_view_bindings);

    auto py_priming_engine_view_bindings = py::class_<gridfire::engine::NetworkPrimingEngineView, gridfire::engine::DefinedEngineView>(m, "NetworkPrimingEngineView");
    py_priming_engine_view_bindings.def(py::init<const std::string&, gridfire::engine::GraphEngine&>(),
        py::arg("primingSymbol"),
        py::arg("baseEngine"),
        "Construct a priming engine view with a priming symbol and a base engine.");
    py_priming_engine_view_bindings.def(py::init<const fourdst::atomic::Species&, gridfire::engine::GraphEngine&>(),
        py::arg("primingSpecies"),
        py::arg("baseEngine"),
        "Construct a priming engine view with a priming species and a base engine.");
    py_priming_engine_view_bindings.def("getBaseEngine", &gridfire::engine::NetworkPrimingEngineView::getBaseEngine,
        "Get the base engine associated with this priming engine view.");

    registerDynamicEngineDefs<gridfire::engine::NetworkPrimingEngineView, gridfire::engine::DefinedEngineView>(py_priming_engine_view_bindings);

    auto py_adaptive_engine_view_bindings = py::class_<gridfire::engine::AdaptiveEngineView, gridfire::engine::DynamicEngine>(m, "AdaptiveEngineView");
    py_adaptive_engine_view_bindings.def(py::init<gridfire::engine::DynamicEngine&>(),
        py::arg("baseEngine"),
        "Construct an adaptive engine view with a base engine.");
    py_adaptive_engine_view_bindings.def("getBaseEngine",
        &gridfire::engine::AdaptiveEngineView::getBaseEngine,
        "Get the base engine associated with this adaptive engine view."
    );

    registerDynamicEngineDefs<gridfire::engine::AdaptiveEngineView, gridfire::engine::DynamicEngine>(py_adaptive_engine_view_bindings);

    auto py_multiscale_engine_view_bindings = py::class_<gridfire::engine::MultiscalePartitioningEngineView, gridfire::engine::DynamicEngine>(m, "MultiscalePartitioningEngineView");
    py_multiscale_engine_view_bindings.def(py::init<gridfire::engine::GraphEngine&>(),
        py::arg("baseEngine"),
        "Construct a multiscale partitioning engine view with a base engine."
    );
    py_multiscale_engine_view_bindings.def("getBaseEngine",
        &gridfire::engine::MultiscalePartitioningEngineView::getBaseEngine,
        "Get the base engine associated with this multiscale partitioning engine view."
    );
    py_multiscale_engine_view_bindings.def("partitionNetwork",
        &gridfire::engine::MultiscalePartitioningEngineView::partitionNetwork,
        py::arg("netIn"),
        "Partition the network based on species timescales and connectivity.");
    py_multiscale_engine_view_bindings.def("partitionNetwork",
        py::overload_cast<const gridfire::NetIn&>(&gridfire::engine::MultiscalePartitioningEngineView::partitionNetwork),
        py::arg("netIn"),
        "Partition the network based on a NetIn object."
    );
    py_multiscale_engine_view_bindings.def("exportToDot",
    &gridfire::engine::MultiscalePartitioningEngineView::exportToDot,
        py::arg("filename"),
        py::arg("comp"),
        py::arg("T9"),
        py::arg("rho"),
        "Export the network to a DOT file for visualization."
    );
    py_multiscale_engine_view_bindings.def("getFastSpecies",
        &gridfire::engine::MultiscalePartitioningEngineView::getFastSpecies,
        "Get the list of fast species in the network."
    );
    py_multiscale_engine_view_bindings.def("getDynamicSpecies",
        &gridfire::engine::MultiscalePartitioningEngineView::getDynamicSpecies,
        "Get the list of dynamic species in the network."
    );
    py_multiscale_engine_view_bindings.def("involvesSpecies",
        &gridfire::engine::MultiscalePartitioningEngineView::involvesSpecies,
        py::arg("species"),
        "Check if a given species is involved in the network (in either the algebraic or dynamic set)."
    );
    py_multiscale_engine_view_bindings.def("involvesSpeciesInQSE",
        &gridfire::engine::MultiscalePartitioningEngineView::involvesSpeciesInQSE,
        py::arg("species"),
        "Check if a given species is involved in the network's algebraic set."
    );
    py_multiscale_engine_view_bindings.def("involvesSpeciesInDynamic",
        &gridfire::engine::MultiscalePartitioningEngineView::involvesSpeciesInDynamic,
        py::arg("species"),
        "Check if a given species is involved in the network's dynamic set."
    );
    py_multiscale_engine_view_bindings.def("getNormalizedEquilibratedComposition",
        &gridfire::engine::MultiscalePartitioningEngineView::getNormalizedEquilibratedComposition,
        py::arg("comp"),
        py::arg("T9"),
        py::arg("rho"),
        "Get the normalized equilibrated composition for the algebraic species."
    );

    registerDynamicEngineDefs<gridfire::engine::MultiscalePartitioningEngineView, gridfire::engine::DynamicEngine>(
        py_multiscale_engine_view_bindings
    );

}









