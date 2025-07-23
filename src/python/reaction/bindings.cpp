#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc if needed directly

#include <string_view>
#include <vector>

#include "bindings.h"

#include "gridfire/reaction/reaction.h"
#include "gridfire/reaction/reaclib.h"

namespace py = pybind11;


void register_reaction_bindings(py::module &m) {
    py::class_<gridfire::reaction::RateCoefficientSet>(m, "RateCoefficientSet")
        .def(py::init<double, double, double, double, double, double, double>(),
            py::arg("a0"), py::arg("a1"), py::arg("a2"), py::arg("a3"),
            py::arg("a4"), py::arg("a5"), py::arg("a6"),
            "Construct a RateCoefficientSet with the given parameters."
        );

    using fourdst::atomic::Species;
    py::class_<gridfire::reaction::Reaction>(m, "Reaction")
        .def(py::init<const std::string_view, const std::string_view, int, const std::vector<Species>&, const std::vector<Species>&, double, std::string_view, gridfire::reaction::RateCoefficientSet, bool>(),
            py::arg("id"), py::arg("peName"), py::arg("chapter"),
            py::arg("reactants"), py::arg("products"), py::arg("qValue"),
            py::arg("label"), py::arg("sets"), py::arg("reverse") = false,
            "Construct a Reaction with the given parameters.")
        .def("calculate_rate", static_cast<double (gridfire::reaction::Reaction::*)(double) const>(&gridfire::reaction::Reaction::calculate_rate),
             py::arg("T9"), "Calculate the reaction rate at a given temperature T9 (in units of 10^9 K).")
        .def("peName", &gridfire::reaction::Reaction::peName,
             "Get the reaction name in (projectile, ejectile) notation (e.g., 'p(p,g)d').")
        .def("chapter", &gridfire::reaction::Reaction::chapter,
            "Get the REACLIB chapter number defining the reaction structure.")
        .def("sourceLabel", &gridfire::reaction::Reaction::sourceLabel,
            "Get the source label for the rate data (e.g., 'wc12w', 'st08').")
        .def("rateCoefficients", &gridfire::reaction::Reaction::rateCoefficients,
            "get the set of rate coefficients.")
        .def("contains", &gridfire::reaction::Reaction::contains,
             py::arg("species"), "Check if the reaction contains a specific species.")
        .def("contains_reactant", &gridfire::reaction::Reaction::contains_reactant,
            "Check if the reaction contains a specific reactant species.")
        .def("contains_product", &gridfire::reaction::Reaction::contains_product,
            "Check if the reaction contains a specific product species.")
        .def("all_species", &gridfire::reaction::Reaction::all_species,
            "Get all species involved in the reaction (both reactants and products) as a set.")
        .def("reactant_species", &gridfire::reaction::Reaction::reactant_species,
            "Get the reactant species of the reaction as a set.")
        .def("product_species", &gridfire::reaction::Reaction::product_species,
            "Get the product species of the reaction as a set.")
        .def("num_species", &gridfire::reaction::Reaction::num_species,
            "Count the number of species in the reaction.")
        .def("stoichiometry", static_cast<int (gridfire::reaction::Reaction::*)(const Species&) const>(&gridfire::reaction::Reaction::stoichiometry),
            py::arg("species"),
            "Get the stoichiometry of the reaction as a map from species to their coefficients.")
        .def("stoichiometry", static_cast<std::unordered_map<Species, int> (gridfire::reaction::Reaction::*)() const>(&gridfire::reaction::Reaction::stoichiometry),
             "Get the stoichiometry of the reaction as a map from species to their coefficients.")
        .def("id", &gridfire::reaction::Reaction::id,
             "Get the unique identifier of the reaction.")
        .def("qValue", &gridfire::reaction::Reaction::qValue,
            "Get the Q-value of the reaction in MeV.")
        .def("reactants", &gridfire::reaction::Reaction::reactants,
            "Get a list of reactant species in the reaction.")
        .def("products", &gridfire::reaction::Reaction::products,
            "Get a list of product species in the reaction.")
        .def("is_reverse", &gridfire::reaction::Reaction::is_reverse,
            "Check if this is a reverse reaction rate.")
        .def("excess_energy", &gridfire::reaction::Reaction::excess_energy,
            "Calculate the excess energy from the mass difference of reactants and products.")
        .def("__eq__", &gridfire::reaction::Reaction::operator==,
            "Equality operator for reactions based on their IDs.")
        .def("__neq__", &gridfire::reaction::Reaction::operator!=,
            "Inequality operator for reactions based on their IDs.")
        .def("hash", &gridfire::reaction::Reaction::hash,
             py::arg("seed") = 0,
             "Compute a hash for the reaction based on its ID.")
        .def("__repr__", [](const gridfire::reaction::Reaction& self) {
            std::stringstream ss;
            ss << self; // Use the existing operator<< for Reaction
            return ss.str();
        });

    py::class_<gridfire::reaction::LogicalReaction, gridfire::reaction::Reaction>(m, "LogicalReaction")
        .def(py::init<const std::vector<gridfire::reaction::Reaction>>(),
             py::arg("reactions"),
             "Construct a LogicalReaction from a vector of Reaction objects.")
        .def("add_reaction", &gridfire::reaction::LogicalReaction::add_reaction,
             py::arg("reaction"),
             "Add another Reaction source to this logical reaction.")
        .def("size", &gridfire::reaction::LogicalReaction::size,
             "Get the number of source rates contributing to this logical reaction.")
        .def("__len__", &gridfire::reaction::LogicalReaction::size,
            "Overload len() to return the number of source rates.")
        .def("sources", &gridfire::reaction::LogicalReaction::sources,
            "Get the list of source labels for the aggregated rates.")
        .def("calculate_rate", static_cast<double (gridfire::reaction::LogicalReaction::*)(double) const>(&gridfire::reaction::LogicalReaction::calculate_rate),
             py::arg("T9"), "Calculate the reaction rate at a given temperature T9 (in units of 10^9 K).")
        .def("calculate_forward_rate_log_derivative", &gridfire::reaction::LogicalReaction::calculate_forward_rate_log_derivative,
             py::arg("T9"), "Calculate the forward rate log derivative at a given temperature T9 (in units of 10^9 K).");

    py::class_<gridfire::reaction::LogicalReactionSet>(m, "LogicalReactionSet")
        .def(py::init<const std::vector<gridfire::reaction::LogicalReaction>>(),
            py::arg("reactions"),
            "Construct a LogicalReactionSet from a vector of LogicalReaction objects.")
        .def(py::init<>(),
            "Default constructor for an empty LogicalReactionSet.")
        .def(py::init<const gridfire::reaction::LogicalReactionSet&>(),
            py::arg("other"),
            "Copy constructor for LogicalReactionSet.")
        .def("add_reaction", &gridfire::reaction::LogicalReactionSet::add_reaction,
            py::arg("reaction"),
            "Add a LogicalReaction to the set.")
        .def("remove_reaction", &gridfire::reaction::LogicalReactionSet::remove_reaction,
            py::arg("reaction"),
            "Remove a LogicalReaction from the set.")
        .def("contains", py::overload_cast<const std::string_view&>(&gridfire::reaction::LogicalReactionSet::contains, py::const_),
            py::arg("id"),
            "Check if the set contains a specific LogicalReaction.")
        .def("contains", py::overload_cast<const gridfire::reaction::Reaction&>(&gridfire::reaction::LogicalReactionSet::contains, py::const_),
            py::arg("reaction"),
            "Check if the set contains a specific Reaction.")
        .def("size", &gridfire::reaction::LogicalReactionSet::size,
            "Get the number of LogicalReactions in the set.")
        .def("__len__", &gridfire::reaction::LogicalReactionSet::size,
            "Overload len() to return the number of LogicalReactions.")
        .def("clear", &gridfire::reaction::LogicalReactionSet::clear,
            "Remove all LogicalReactions from the set.")
        .def("containes_species", &gridfire::reaction::LogicalReactionSet::contains_species,
            py::arg("species"),
            "Check if any reaction in the set involves the given species.")
        .def("contains_reactant", &gridfire::reaction::LogicalReactionSet::contains_reactant,
            py::arg("species"),
            "Check if any reaction in the set has the species as a reactant.")
        .def("contains_product", &gridfire::reaction::LogicalReactionSet::contains_product,
            py::arg("species"),
            "Check if any reaction in the set has the species as a product.")
        .def("__getitem__", py::overload_cast<size_t>(&gridfire::reaction::LogicalReactionSet::operator[], py::const_),
            py::arg("index"),
            "Get a LogicalReaction by index.")
        .def("__getitem___", py::overload_cast<const std::string_view&>(&gridfire::reaction::LogicalReactionSet::operator[], py::const_),
            py::arg("id"),
            "Get a LogicalReaction by its ID.")
        .def("__eq__", &gridfire::reaction::LogicalReactionSet::operator==,
            py::arg("LogicalReactionSet"),
            "Equality operator for LogicalReactionSets based on their contents.")
        .def("__ne__", &gridfire::reaction::LogicalReactionSet::operator!=,
            py::arg("LogicalReactionSet"),
            "Inequality operator for LogicalReactionSets based on their contents.")
        .def("hash", &gridfire::reaction::LogicalReactionSet::hash,
             py::arg("seed") = 0,
             "Compute a hash for the LogicalReactionSet based on its contents."
        )
        .def("__repr__", [](const gridfire::reaction::LogicalReactionSet& self) {
            std::stringstream ss;
            ss << self;
            return ss.str();
        })
        .def("getReactionSetSpecies", &gridfire::reaction::LogicalReactionSet::getReactionSetSpecies,
             "Get all species involved in the reactions of the set as a set of Species objects.");

        m.def("packReactionSetToLogicalReactionSet",
            &gridfire::reaction::packReactionSetToLogicalReactionSet,
            py::arg("reactionSet"),
            "Convert a ReactionSet to a LogicalReactionSet by aggregating reactions with the same peName."
        );

        m.def("get_all_reactions", &gridfire::reaclib::get_all_reactions,
            "Get all reactions from the REACLIB database.");
}
