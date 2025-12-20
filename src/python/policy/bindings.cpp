#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <memory>
#include "bindings.h"
#include "trampoline/py_policy.h"

#include "gridfire/policy/policy.h"


namespace py = pybind11;

namespace {

    template <typename T>
    concept IsReactionChainPolicy = std::is_base_of_v<gridfire::policy::ReactionChainPolicy, T>;

    template <typename T>
    concept IsNetworkPolicy = std::is_base_of_v<gridfire::policy::NetworkPolicy, T>;

    template<IsReactionChainPolicy T, IsReactionChainPolicy BaseT>
    void registerReactionChainPolicyDefs(py::class_<T, BaseT>& pyClass) {
        pyClass.def(
            "get_reactions",
            &T::get_reactions,
            "Get the ReactionSet representing this reaction chain."
        )
        .def(
            "contains",
            py::overload_cast<const std::string&>(&T::contains, py::const_),
            py::arg("id"),
            "Check if the reaction chain contains a reaction with the given ID."
        )
        .def(
            "contains",
            py::overload_cast<const gridfire::reaction::Reaction&>(&T::contains, py::const_),
            py::arg("reaction"),
            "Check if the reaction chain contains the given reaction."
        )
        .def(
            "name",
            &T::name,
            "Get the name of the reaction chain policy."
        )
        .def(
            "hash",
            &T::hash,
            py::arg("seed"),
            "Compute a hash value for the reaction chain policy."
        )
        .def(
            "__eq__",
            &T::operator==,
            py::arg("other"),
            "Check equality with another ReactionChainPolicy."
        )
        .def(
            "__ne__",
            &T::operator!=,
            py::arg("other"),
            "Check inequality with another ReactionChainPolicy."
        )
        .def("__hash__", [](const T &self) {
            return self.hash(0);
            }
        )
        .def("__repr__", [](const T &self) {
            std::stringstream ss;
            ss << self;
            return ss.str();
        });
    }

    template<IsNetworkPolicy T, IsNetworkPolicy BaseT>
    void registerNetworkPolicyDefs(py::class_<T, BaseT> pyClass) {
        pyClass.def(
            "name",
            &T::name,
            "Get the name of the network policy."
        )
        .def(
            "get_seed_species",
            &T::get_seed_species,
            "Get the set of seed species required by the network policy."
        )
        .def(
            "get_seed_reactions",
            &T::get_seed_reactions,
            "Get the set of seed reactions required by the network policy."
        )
        .def(
            "get_status",
            &T::get_status,
            "Get the current status of the network policy."
        )
        .def(
            "get_engine_types_stack",
            &T::get_engine_types_stack,
            "Get the types of engines in the stack constructed by the network policy."
        )
        .def(
            "construct",
            &T::construct,
            "Construct the network according to the policy."

        )
        .def(
            "get_engine_stack",
            [](const T &self) {
                const auto& stack = self.get_engine_stack();
                std::vector<gridfire::engine::DynamicEngine*> engine_ptrs;
                engine_ptrs.reserve(stack.size());
                for (const auto& engine_uptr : stack) {
                    engine_ptrs.push_back(engine_uptr.get());
                }

                return engine_ptrs;
            },
            py::return_value_policy::reference_internal
        )
        .def(
            "get_stack_scratch_blob",
            &T::get_stack_scratch_blob
        )
        .def(
            "get_partition_function",
            &T::get_partition_function
        );
    }
}


void register_policy_bindings(pybind11::module &m) {
    register_reaction_chain_policy_bindings(m);
    register_network_policy_bindings(m);
}

void register_reaction_chain_policy_bindings(pybind11::module &m) {
    using namespace gridfire::policy;

    py::class_<ReactionChainPolicy, PyReactionChainPolicy> py_reactionChainPolicy(m, "ReactionChainPolicy");
    py::class_<MultiReactionChainPolicy, ReactionChainPolicy> py_multiChainPolicy(m, "MultiReactionChainPolicy");
    py::class_<TemperatureDependentChainPolicy, ReactionChainPolicy> py_tempDepChainPolicy(m, "TemperatureDependentChainPolicy");


    py::class_<ProtonProtonIChainPolicy, TemperatureDependentChainPolicy> py_ppI(m, "ProtonProtonIChainPolicy");
    py_ppI.def(py::init<>());
    py_ppI.def("name", &ProtonProtonIChainPolicy::name, "Get the name of the reaction chain policy.");

    py::class_<ProtonProtonIIChainPolicy, TemperatureDependentChainPolicy> py_ppII(m, "ProtonProtonIIChainPolicy");
    py_ppII.def(py::init<>());
    py_ppII.def("name", &ProtonProtonIIChainPolicy::name, "Get the name of the reaction chain policy.");

    py::class_<ProtonProtonIIIChainPolicy, TemperatureDependentChainPolicy> py_ppIII(m, "ProtonProtonIIIChainPolicy");
    py_ppIII.def(py::init<>());
    py_ppIII.def("name", &ProtonProtonIIIChainPolicy::name, "Get the name of the reaction chain policy.");

    py::class_<ProtonProtonChainPolicy, MultiReactionChainPolicy> py_ppChain(m, "ProtonProtonChainPolicy");
    py_ppChain.def(py::init<>());
    py_ppChain.def("name", &ProtonProtonChainPolicy::name, "Get the name of the reaction chain policy.");

    registerReactionChainPolicyDefs(py_ppI);
    registerReactionChainPolicyDefs(py_ppII);
    registerReactionChainPolicyDefs(py_ppIII);
    registerReactionChainPolicyDefs(py_ppChain);

    py::class_<CNOIChainPolicy, TemperatureDependentChainPolicy> py_cnoI(m, "CNOIChainPolicy");
    py_cnoI.def(py::init<>());
    py_cnoI.def("name", &CNOIChainPolicy::name, "Get the name of the reaction chain policy.");

    py::class_<CNOIIChainPolicy, TemperatureDependentChainPolicy> py_cnoII(m, "CNOIIChainPolicy");
    py_cnoII.def(py::init<>());
    py_cnoII.def("name", &CNOIIChainPolicy::name, "Get the name of the reaction chain policy.");

    py::class_<CNOIIIChainPolicy, TemperatureDependentChainPolicy> py_cnoIII(m, "CNOIIIChainPolicy");
    py_cnoIII.def(py::init<>());
    py_cnoIII.def("name", &CNOIIIChainPolicy::name, "Get the name of the reaction chain policy.");

    py::class_<CNOIVChainPolicy, TemperatureDependentChainPolicy> py_cnoIV(m, "CNOIVChainPolicy");
    py_cnoIV.def(py::init<>());
    py_cnoIV.def("name", &CNOIVChainPolicy::name, "Get the name of the reaction chain policy.");

    py::class_<CNOChainPolicy, MultiReactionChainPolicy> py_cnoChain(m, "CNOChainPolicy");
    py_cnoChain.def(py::init<>());
    py_cnoChain.def("name", &CNOChainPolicy::name, "Get the name of the reaction chain policy.");

    registerReactionChainPolicyDefs(py_cnoI);
    registerReactionChainPolicyDefs(py_cnoII);
    registerReactionChainPolicyDefs(py_cnoIII);
    registerReactionChainPolicyDefs(py_cnoIV);
    registerReactionChainPolicyDefs(py_cnoChain);

    py::class_<HotCNOIChainPolicy, TemperatureDependentChainPolicy> py_hotCNOI(m, "HotCNOIChainPolicy");
    py_hotCNOI.def(py::init<>());
    py_hotCNOI.def("name", &HotCNOIChainPolicy::name, "Get the name of the reaction chain policy.");

    py::class_<HotCNOIIChainPolicy, TemperatureDependentChainPolicy> py_hotCNOII(m, "HotCNOIIChainPolicy");
    py_hotCNOII.def(py::init<>());
    py_hotCNOII.def("name", &HotCNOIIChainPolicy::name, "Get the name of the reaction chain policy.");

    py::class_<HotCNOIIIChainPolicy, TemperatureDependentChainPolicy> py_hotCNOIII(m, "HotCNOIIIChainPolicy");
    py_hotCNOIII.def(py::init<>());
    py_hotCNOIII.def("name", &HotCNOIIIChainPolicy::name, "Get the name of the reaction chain policy.");

    py::class_<HotCNOChainPolicy, MultiReactionChainPolicy> py_hotCNOChain(m, "HotCNOChainPolicy");
    py_hotCNOChain.def(py::init<>());
    py_hotCNOChain.def("name", &HotCNOChainPolicy::name, "Get the name of the reaction chain policy.");

    registerReactionChainPolicyDefs(py_hotCNOI);
    registerReactionChainPolicyDefs(py_hotCNOII);
    registerReactionChainPolicyDefs(py_hotCNOIII);
    registerReactionChainPolicyDefs(py_hotCNOChain);

    py::class_<TripleAlphaChainPolicy, TemperatureDependentChainPolicy> py_tripleAlpha(m, "TripleAlphaChainPolicy");
    py_tripleAlpha.def(py::init<>());
    py_tripleAlpha.def("name", &TripleAlphaChainPolicy::name, "Get the name of the reaction chain policy.");

    registerReactionChainPolicyDefs(py_tripleAlpha);

    py::class_<MainSequenceReactionChainPolicy, MultiReactionChainPolicy> py_mainSeq(m, "MainSequenceReactionChainPolicy");
    py_mainSeq.def(py::init<>());
    py_mainSeq.def("name", &MainSequenceReactionChainPolicy::name, "Get the name of the reaction chain policy.");

    registerReactionChainPolicyDefs(py_mainSeq);

}

void register_network_policy_bindings(pybind11::module &m) {
    py::enum_<gridfire::policy::NetworkPolicyStatus>(m, "NetworkPolicyStatus")
        .value("UNINITIALIZED", gridfire::policy::NetworkPolicyStatus::UNINITIALIZED)
        .value("INITIALIZED_UNVERIFIED", gridfire::policy::NetworkPolicyStatus::INITIALIZED_UNVERIFIED)
        .value("MISSING_KEY_REACTION", gridfire::policy::NetworkPolicyStatus::MISSING_KEY_REACTION)
        .value("MISSING_KEY_SPECIES", gridfire::policy::NetworkPolicyStatus::MISSING_KEY_SPECIES)
        .value("INITIALIZED_VERIFIED", gridfire::policy::NetworkPolicyStatus::INITIALIZED_VERIFIED)
        .export_values();

    m.def("network_policy_status_to_string",
        &gridfire::policy::NetworkPolicyStatusToString,
        py::arg("status"),
        "Convert a NetworkPolicyStatus enum value to its string representation."
    );

    py::class_<gridfire::policy::ConstructionResults>(m, "ConstructionResults")
        .def_property_readonly("engine",
            [](const gridfire::policy::ConstructionResults &self) -> const gridfire::engine::DynamicEngine& {
                return self.engine;
            },
            py::return_value_policy::reference
        )
        .def_property_readonly("scratch_blob",
            [](const gridfire::policy::ConstructionResults &self) {
                return self.scratch_blob.get();
            },
            py::return_value_policy::reference_internal
        );

    py::class_<gridfire::policy::NetworkPolicy, PyNetworkPolicy> py_networkPolicy(m, "NetworkPolicy");
    py::class_<gridfire::policy::MainSequencePolicy, gridfire::policy::NetworkPolicy> py_mainSeqPolicy(m, "MainSequencePolicy");
    py_mainSeqPolicy.def(
        py::init<const fourdst::composition::Composition&>(),
        py::arg("composition"),
        "Construct MainSequencePolicy from an existing composition."
    );
    py_mainSeqPolicy.def(
        py::init<std::vector<fourdst::atomic::Species>, const std::vector<double>&>(),
        py::arg("seed_species"),
        py::arg("mass_fractions"),
        "Construct MainSequencePolicy from seed species and mass fractions."
    );

    registerNetworkPolicyDefs(py_mainSeqPolicy);
}
