#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc if needed directly

#include <iostream>
#include <memory>

#include "bindings.h"

#include "gridfire/partition/partition.h"

PYBIND11_DECLARE_HOLDER_TYPE(T, std::unique_ptr<T>, true) // Declare unique_ptr as a holder type for pybind11

#include "trampoline/py_partition.h"


namespace py = pybind11;


void register_partition_bindings(pybind11::module &m) {
    using PF = gridfire::partition::PartitionFunction;
    py::class_<PF, PyPartitionFunction>(m, "PartitionFunction");

    register_partition_types_bindings(m);
    register_ground_state_partition_bindings(m);
    register_rauscher_thielemann_partition_data_record_bindings(m);
    register_rauscher_thielemann_partition_bindings(m);

    register_composite_partition_bindings(m);
}

void register_partition_types_bindings(pybind11::module &m) {
    py::enum_<gridfire::partition::BasePartitionType>(m, "BasePartitionType")
        .value("RauscherThielemann", gridfire::partition::BasePartitionType::RauscherThielemann)
        .value("GroundState", gridfire::partition::BasePartitionType::GroundState)
        .export_values();

    m.def("basePartitionTypeToString", [](gridfire::partition::BasePartitionType type) {
        return gridfire::partition::basePartitionTypeToString[type];
    }, py::arg("type"), "Convert BasePartitionType to string.");

    m.def("stringToBasePartitionType", [](const std::string &typeStr) {
        return gridfire::partition::stringToBasePartitionType[typeStr];
    }, py::arg("typeStr"), "Convert string to BasePartitionType.");
}

void register_ground_state_partition_bindings(pybind11::module &m) {
    using GSPF = gridfire::partition::GroundStatePartitionFunction;
    using PF = gridfire::partition::PartitionFunction;
    py::class_<GSPF, PF>(m, "GroundStatePartitionFunction")
        .def(py::init<>())
        .def("evaluate", &gridfire::partition::GroundStatePartitionFunction::evaluate,
             py::arg("z"), py::arg("a"), py::arg("T9"),
             "Evaluate the ground state partition function for given Z, A, and T9.")
        .def("evaluateDerivative", &gridfire::partition::GroundStatePartitionFunction::evaluateDerivative,
             py::arg("z"), py::arg("a"), py::arg("T9"),
             "Evaluate the derivative of the ground state partition function for given Z, A, and T9.")
        .def("supports", &gridfire::partition::GroundStatePartitionFunction::supports,
             py::arg("z"), py::arg("a"),
             "Check if the ground state partition function supports given Z and A.")
        .def("get_type", &gridfire::partition::GroundStatePartitionFunction::type,
             "Get the type of the partition function (should return 'GroundState').");
}

void register_rauscher_thielemann_partition_data_record_bindings(pybind11::module &m) {
    py::class_<gridfire::partition::record::RauscherThielemannPartitionDataRecord>(m, "RauscherThielemannPartitionDataRecord")
        .def_readonly("z", &gridfire::partition::record::RauscherThielemannPartitionDataRecord::z, "Atomic number")
        .def_readonly("a", &gridfire::partition::record::RauscherThielemannPartitionDataRecord::a, "Mass number")
        .def_readonly("ground_state_spin", &gridfire::partition::record::RauscherThielemannPartitionDataRecord::ground_state_spin, "Ground state spin")
        .def_readonly("normalized_g_values", &gridfire::partition::record::RauscherThielemannPartitionDataRecord::normalized_g_values, "Normalized g-values for the first 24 energy levels");
}


void register_rauscher_thielemann_partition_bindings(pybind11::module &m) {
    using RTPF = gridfire::partition::RauscherThielemannPartitionFunction;
    using PF = gridfire::partition::PartitionFunction;
    py::class_<RTPF, PF>(m, "RauscherThielemannPartitionFunction")
        .def(py::init<>())
        .def("evaluate", &gridfire::partition::RauscherThielemannPartitionFunction::evaluate,
             py::arg("z"), py::arg("a"), py::arg("T9"),
             "Evaluate the Rauscher-Thielemann partition function for given Z, A, and T9.")
        .def("evaluateDerivative", &gridfire::partition::RauscherThielemannPartitionFunction::evaluateDerivative,
             py::arg("z"), py::arg("a"), py::arg("T9"),
             "Evaluate the derivative of the Rauscher-Thielemann partition function for given Z, A, and T9.")
        .def("supports", &gridfire::partition::RauscherThielemannPartitionFunction::supports,
             py::arg("z"), py::arg("a"),
             "Check if the Rauscher-Thielemann partition function supports given Z and A.")
        .def("get_type", &gridfire::partition::RauscherThielemannPartitionFunction::type,
             "Get the type of the partition function (should return 'RauscherThielemann').");
}

void register_composite_partition_bindings(pybind11::module &m) {
    py::class_<gridfire::partition::CompositePartitionFunction>(m, "CompositePartitionFunction")
        .def(py::init<const std::vector<gridfire::partition::BasePartitionType>&>(),
             py::arg("partitionFunctions"),
             "Create a composite partition function from a list of base partition types.")
        .def(py::init<const gridfire::partition::CompositePartitionFunction&>(),
             "Copy constructor for CompositePartitionFunction.")
        .def("evaluate", &gridfire::partition::CompositePartitionFunction::evaluate,
             py::arg("z"), py::arg("a"), py::arg("T9"),
             "Evaluate the composite partition function for given Z, A, and T9.")
        .def("evaluateDerivative", &gridfire::partition::CompositePartitionFunction::evaluateDerivative,
             py::arg("z"), py::arg("a"), py::arg("T9"),
             "Evaluate the derivative of the composite partition function for given Z, A, and T9.")
        .def("supports", &gridfire::partition::CompositePartitionFunction::supports,
             py::arg("z"), py::arg("a"),
             "Check if the composite partition function supports given Z and A.")
        .def("get_type", &gridfire::partition::CompositePartitionFunction::type,
             "Get the type of the partition function (should return 'Composite').");
}



