#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly

#include "gridfire/engine/scratchpads/scratchpads.h"

#include "bindings.h"

namespace py = pybind11;
namespace sp = gridfire::engine::scratch;

template<typename... ScratchPadTypes>
void build_state_getter(py::module& m) {

}

void register_scratchpad_types_bindings(pybind11::module &m) {
    py::enum_<sp::ScratchPadType>(m, "ScratchPadType")
        .value("GRAPH_ENGINE_SCRATCHPAD", sp::ScratchPadType::GRAPH_ENGINE_SCRATCHPAD)
        .value("MULTISCALE_PARTITIONING_ENGINE_VIEW_SCRATCHPAD", sp::ScratchPadType::MULTISCALE_PARTITIONING_ENGINE_VIEW_SCRATCHPAD)
        .value("ADAPTIVE_ENGINE_VIEW_SCRATCHPAD", sp::ScratchPadType::ADAPTIVE_ENGINE_VIEW_SCRATCHPAD)
        .value("DEFINED_ENGINE_VIEW_SCRATCHPAD", sp::ScratchPadType::DEFINED_ENGINE_VIEW_SCRATCHPAD)
        .export_values();
}

void register_scratchpad_bindings(pybind11::module_ &m) {
    py::enum_<sp::GraphEngineScratchPad::ADFunRegistrationResult>(m, "ADFunRegistrationResult")
        .value("SUCCESS", sp::GraphEngineScratchPad::ADFunRegistrationResult::SUCCESS)
        .value("ALREADY_REGISTERED", sp::GraphEngineScratchPad::ADFunRegistrationResult::ALREADY_REGISTERED)
        .export_values();

    py::class_<sp::GraphEngineScratchPad>(m, "GraphEngineScratchPad")
        .def(py::init<>())
        .def("initialize", &sp::GraphEngineScratchPad::initialize, py::arg("engine"))
        .def("clone", &sp::GraphEngineScratchPad::clone)
        .def("is_initialized", &sp::GraphEngineScratchPad::is_initialized)
        .def_readonly("most_recent_rhs_calculation", &sp::GraphEngineScratchPad::most_recent_rhs_calculation)
        .def_readonly("local_abundance_cache", &sp::GraphEngineScratchPad::local_abundance_cache)
        .def_readonly("has_initialized", &sp::GraphEngineScratchPad::has_initialized)
        .def_readonly("stepDerivativesCache", &sp::GraphEngineScratchPad::stepDerivativesCache)
        .def_readonly_static("ID", &sp::GraphEngineScratchPad::ID)
        .def("__repr__", [](const sp::GraphEngineScratchPad &self) {
            return std::format("{}", self);
        });

    py::class_<sp::MultiscalePartitioningEngineViewScratchPad>(m, "MultiscalePartitioningEngineViewScratchPad")
        .def(py::init<>())
        .def("initialize", &sp::MultiscalePartitioningEngineViewScratchPad::initialize)
        .def("clone", &sp::MultiscalePartitioningEngineViewScratchPad::clone)
        .def("is_initialized", &sp::MultiscalePartitioningEngineViewScratchPad::is_initialized)
        .def_readonly("qse_groups", &sp::MultiscalePartitioningEngineViewScratchPad::qse_groups)
        .def_readonly("dynamic_species", &sp::MultiscalePartitioningEngineViewScratchPad::dynamic_species)
        .def_readonly("algebraic_species", &sp::MultiscalePartitioningEngineViewScratchPad::algebraic_species)
        .def_readonly("composition_cache", &sp::MultiscalePartitioningEngineViewScratchPad::composition_cache)
        .def_readonly("has_initialized", &sp::MultiscalePartitioningEngineViewScratchPad::has_initialized)
        .def_readonly_static("ID", &sp::MultiscalePartitioningEngineViewScratchPad::ID)
        .def("__repr__", [](const sp::MultiscalePartitioningEngineViewScratchPad &self) {
            return std::format("{}", self);
        });

    py::class_<sp::AdaptiveEngineViewScratchPad>(m, "AdaptiveEngineViewScratchPad")
        .def(py::init<>())
        .def("initialize", &sp::AdaptiveEngineViewScratchPad::initialize)
        .def("clone", &sp::AdaptiveEngineViewScratchPad::clone)
        .def("is_initialized", &sp::AdaptiveEngineViewScratchPad::is_initialized)
        .def_readonly("active_species", &sp::AdaptiveEngineViewScratchPad::active_species)
        .def_readonly("active_reactions", &sp::AdaptiveEngineViewScratchPad::active_reactions)
        .def_readonly("has_initialized", &sp::AdaptiveEngineViewScratchPad::has_initialized)
        .def_readonly_static("ID", &sp::AdaptiveEngineViewScratchPad::ID)
        .def("__repr__", [](const sp::AdaptiveEngineViewScratchPad &self) {
            return std::format("{}", self);
        });

    py::class_<sp::DefinedEngineViewScratchPad>(m, "DefinedEngineViewScratchPad")
        .def(py::init<>())
        .def("clone", &sp::DefinedEngineViewScratchPad::clone)
        .def("is_initialized", &sp::DefinedEngineViewScratchPad::is_initialized)
        .def_readonly("active_species", &sp::DefinedEngineViewScratchPad::active_species)
        .def_readonly("active_reactions", &sp::DefinedEngineViewScratchPad::active_reactions)
        .def_readonly("species_index_map", &sp::DefinedEngineViewScratchPad::species_index_map)
        .def_readonly("reaction_index_map", &sp::DefinedEngineViewScratchPad::reaction_index_map)
        .def_readonly("has_initialized", &sp::DefinedEngineViewScratchPad::has_initialized)
        .def_readonly_static("ID", &sp::DefinedEngineViewScratchPad::ID)
        .def("__repr__", [](const sp::DefinedEngineViewScratchPad &self) {
            return std::format("{}", self);
        });
}

void register_state_blob_bindings(pybind11::module_ &m) {
    py::enum_<sp::StateBlob::Error>(m, "StateBlobError")
        .value("SCRATCHPAD_OUT_OF_BOUNDS", sp::StateBlob::Error::SCRATCHPAD_OUT_OF_BOUNDS)
        .value("SCRATCHPAD_NOT_FOUND", sp::StateBlob::Error::SCRATCHPAD_NOT_FOUND)
        .value("SCRATCHPAD_BAD_CAST", sp::StateBlob::Error::SCRATCHPAD_BAD_CAST)
        .value("SCRATCHPAD_NOT_INITIALIZED", sp::StateBlob::Error::SCRATCHPAD_NOT_INITIALIZED)
        .value("SCRATCHPAD_TYPE_COLLISION", sp::StateBlob::Error::SCRATCHPAD_TYPE_COLLISION)
        .value("SCRATCHPAD_UNKNOWN_ERROR", sp::StateBlob::Error::SCRATCHPAD_UNKNOWN_ERROR)
        .export_values();

    py::class_<sp::StateBlob>(m, "StateBlob")
        .def(py::init<>())
        .def("enroll", [](sp::StateBlob &self, const sp::ScratchPadType type) {
            switch (type) {
                case sp::ScratchPadType::GRAPH_ENGINE_SCRATCHPAD:
                    self.enroll<sp::GraphEngineScratchPad>();
                    break;
                case sp::ScratchPadType::MULTISCALE_PARTITIONING_ENGINE_VIEW_SCRATCHPAD:
                    self.enroll<sp::MultiscalePartitioningEngineViewScratchPad>();
                    break;
                case sp::ScratchPadType::ADAPTIVE_ENGINE_VIEW_SCRATCHPAD:
                    self.enroll<sp::AdaptiveEngineViewScratchPad>();
                    break;
                case sp::ScratchPadType::DEFINED_ENGINE_VIEW_SCRATCHPAD:
                    self.enroll<sp::DefinedEngineViewScratchPad>();
                    break;
                default:
                    throw std::invalid_argument("Unknown ScratchPadType for enrollment.");
            }
        })
    .def("get", [](const sp::StateBlob &self, const sp::ScratchPadType type) {
            auto result = self.get(type);
            if (!result.has_value()) {
                throw std::runtime_error("Error retrieving scratchpad: " + sp::StateBlob::error_to_string(result.error()));
            }
            return result.value();
        },
        pybind11::return_value_policy::reference_internal
    )
    .def("clone_structure", &sp::StateBlob::clone_structure)
    .def("get_registered_scratchpads", &sp::StateBlob::get_registered_scratchpads)
    .def("get_status", [](const sp::StateBlob &self, const sp::ScratchPadType type) -> sp::StateBlob::ScratchPadStatus {
            switch (type) {
                case sp::ScratchPadType::GRAPH_ENGINE_SCRATCHPAD:
                    return self.get_status<sp::GraphEngineScratchPad>();
                case sp::ScratchPadType::MULTISCALE_PARTITIONING_ENGINE_VIEW_SCRATCHPAD:
                    return self.get_status<sp::MultiscalePartitioningEngineViewScratchPad>();
                case sp::ScratchPadType::ADAPTIVE_ENGINE_VIEW_SCRATCHPAD:
                    return self.get_status<sp::AdaptiveEngineViewScratchPad>();
                case sp::ScratchPadType::DEFINED_ENGINE_VIEW_SCRATCHPAD:
                    return self.get_status<sp::DefinedEngineViewScratchPad>();
                default:
                    throw std::invalid_argument("Unknown ScratchPadType for status retrieval.");
            }
    })
    .def("get_status_map", &sp::StateBlob::get_status_map)
    .def_static("error_to_string", &sp::StateBlob::error_to_string)
    .def("__repr__", [](const sp::StateBlob &self) {
        return std::format("{}", self);
    });
}



