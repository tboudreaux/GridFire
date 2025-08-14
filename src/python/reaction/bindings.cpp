#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly

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
    py::class_<gridfire::reaction::ReaclibReaction>(m, "ReaclibReaction")
        .def(
            py::init<
                const std::string_view,
                const std::string_view,
                int,
                const std::vector<Species>&,
                const std::vector<Species>&,
                double, std::string_view,
                gridfire::reaction::RateCoefficientSet,
                bool
            >(),
            py::arg("id"),
            py::arg("peName"),
            py::arg("chapter"),
            py::arg("reactants"),
            py::arg("products"),
            py::arg("qValue"),
            py::arg("label"),
            py::arg("sets"),
            py::arg("reverse") = false,
            "Construct a Reaction with the given parameters."
        )
        .def(
            "calculate_rate",
             [](const gridfire::reaction::ReaclibReaction& self, const double T9, const double rho, const std::vector<double>& Y) -> double {
                 return self.calculate_rate(T9, rho, Y);
             },
             py::arg("T9"),
             py::arg("rho"),
             py::arg("Y"),
             "Calculate the reaction rate at a given temperature T9 (in units of 10^9 K)."
        )
        .def(
            "peName",
            &gridfire::reaction::ReaclibReaction::peName,
             "Get the reaction name in (projectile, ejectile) notation (e.g., 'p(p,g)d')."
        )
        .def(
            "chapter",
            &gridfire::reaction::ReaclibReaction::chapter,
            "Get the REACLIB chapter number defining the reaction structure."
        )
        .def(
            "sourceLabel",
            &gridfire::reaction::ReaclibReaction::sourceLabel,
            "Get the source label for the rate data (e.g., 'wc12w', 'st08')."
        )
        .def(
            "rateCoefficients",
            &gridfire::reaction::ReaclibReaction::rateCoefficients,
            "get the set of rate coefficients."
        )
        .def(
            "contains",
            &gridfire::reaction::ReaclibReaction::contains,
             py::arg("species"),
             "Check if the reaction contains a specific species."
        )
        .def(
            "contains_reactant",
            &gridfire::reaction::ReaclibReaction::contains_reactant,
            "Check if the reaction contains a specific reactant species."
        )
        .def(
            "contains_product",
            &gridfire::reaction::ReaclibReaction::contains_product,
            "Check if the reaction contains a specific product species."
        )
        .def(
            "all_species",
            &gridfire::reaction::ReaclibReaction::all_species,
            "Get all species involved in the reaction (both reactants and products) as a set."
        )
        .def(
            "reactant_species",
            &gridfire::reaction::ReaclibReaction::reactant_species,
            "Get the reactant species of the reaction as a set."
        )
        .def(
            "product_species",
            &gridfire::reaction::ReaclibReaction::product_species,
            "Get the product species of the reaction as a set."
        )
        .def(
            "num_species",
            &gridfire::reaction::ReaclibReaction::num_species,
            "Count the number of species in the reaction."
        )
        .def(
            "stoichiometry",
            [](const gridfire::reaction::ReaclibReaction& self, const Species& species) -> int {
                return self.stoichiometry(species);
            },
            py::arg("species"),
            "Get the stoichiometry of the reaction as a map from species to their coefficients."
        )
        .def(
            "stoichiometry",
            [](const gridfire::reaction::ReaclibReaction& self) -> std::unordered_map<Species, int> {
                return self.stoichiometry();
            },
             "Get the stoichiometry of the reaction as a map from species to their coefficients."
        )
        .def(
            "id",
            &gridfire::reaction::ReaclibReaction::id,
             "Get the unique identifier of the reaction."
        )
        .def(
            "qValue",
            &gridfire::reaction::ReaclibReaction::qValue,
            "Get the Q-value of the reaction in MeV."
        )
        .def(
            "reactants",
            &gridfire::reaction::ReaclibReaction::reactants,
            "Get a list of reactant species in the reaction."
        )
        .def(
            "products",
            &gridfire::reaction::ReaclibReaction::products,
            "Get a list of product species in the reaction."
        )
        .def(
            "is_reverse",
            &gridfire::reaction::ReaclibReaction::is_reverse,
            "Check if this is a reverse reaction rate."
        )
        .def(
            "excess_energy",
            &gridfire::reaction::ReaclibReaction::excess_energy,
            "Calculate the excess energy from the mass difference of reactants and products."
        )
        .def(
            "__eq__",
            &gridfire::reaction::ReaclibReaction::operator==,
            "Equality operator for reactions based on their IDs."
        )
        .def(
            "__neq__",
            &gridfire::reaction::ReaclibReaction::operator!=,
            "Inequality operator for reactions based on their IDs."
        )
        .def(
            "hash",
            &gridfire::reaction::ReaclibReaction::hash,
             py::arg("seed") = 0,
             "Compute a hash for the reaction based on its ID."
        )
        .def(
            "__repr__",
            [](const gridfire::reaction::ReaclibReaction& self) {
                std::stringstream ss;
                ss << self; // Use the existing operator<< for Reaction
                return ss.str();
            }
        );

    py::class_<gridfire::reaction::LogicalReaclibReaction, gridfire::reaction::ReaclibReaction>(m, "LogicalReaclibReaction")
        .def(
            py::init<const std::vector<gridfire::reaction::Reaction>>(),
             py::arg("reactions"),
             "Construct a LogicalReaclibReaction from a vector of Reaction objects."
        )
        .def(
            "add_reaction",
            &gridfire::reaction::LogicalReaclibReaction::add_reaction,
             py::arg("reaction"),
             "Add another Reaction source to this logical reaction."
        )
        .def(
            "size",
            &gridfire::reaction::LogicalReaclibReaction::size,
             "Get the number of source rates contributing to this logical reaction."
        )
        .def(
            "__len__",
            &gridfire::reaction::LogicalReaclibReaction::size,
            "Overload len() to return the number of source rates."
        )
        .def(
            "sources",
            &gridfire::reaction::LogicalReaclibReaction::sources,
            "Get the list of source labels for the aggregated rates."
        )
        .def(
            "calculate_rate",
            [](const gridfire::reaction::LogicalReaclibReaction& self, const double T9, const double rho, const std::vector<double>& Y) -> double {
                return self.calculate_rate(T9, rho, Y);
            },
            py::arg("T9"),
            "Calculate the reaction rate at a given temperature T9 (in units of 10^9 K)."
        )
        .def(
            "calculate_forward_rate_log_derivative",
            &gridfire::reaction::LogicalReaclibReaction::calculate_forward_rate_log_derivative,
            py::arg("T9"),
            "Calculate the forward rate log derivative at a given temperature T9 (in units of 10^9 K)."
        );

    py::class_<gridfire::reaction::ReactionSet>(m, "ReactionSet")
        // TODO: Fix the constructor to accept a vector of unique ptrs to Reaclib Reactions
        .def(
            py::init<const std::vector<gridfire::reaction::Reaction>>(),
            py::arg("reactions"),
            "Construct a LogicalReactionSet from a vector of LogicalReaclibReaction objects."
        )
        .def(
            py::init<>(),
            "Default constructor for an empty LogicalReactionSet."
        )
        .def(
            py::init<const gridfire::reaction::ReactionSet&>(),
            py::arg("other"),
            "Copy constructor for LogicalReactionSet."
        )
        .def(
            "add_reaction",
            py::overload_cast<const gridfire::reaction::Reaction&>(&gridfire::reaction::ReactionSet::add_reaction),
            py::arg("reaction"),
            "Add a LogicalReaclibReaction to the set."
        )
        .def(
            "remove_reaction",
            &gridfire::reaction::ReactionSet::remove_reaction,
            py::arg("reaction"),
            "Remove a LogicalReaclibReaction from the set."
        )
        .def(
            "contains",
            py::overload_cast<const std::string_view&>(&gridfire::reaction::ReactionSet::contains, py::const_),
            py::arg("id"),
            "Check if the set contains a specific LogicalReaclibReaction."
        )
        .def(
            "contains",
            py::overload_cast<const gridfire::reaction::Reaction&>(&gridfire::reaction::ReactionSet::contains, py::const_),
            py::arg("reaction"),
            "Check if the set contains a specific Reaction."
        )
        .def(
            "size",
            &gridfire::reaction::ReactionSet::size,
            "Get the number of LogicalReactions in the set."
        )
        .def(
            "__len__", &gridfire::reaction::ReactionSet::size,
            "Overload len() to return the number of LogicalReactions."
        )
        .def(
            "clear",
            &gridfire::reaction::ReactionSet::clear,
            "Remove all LogicalReactions from the set."
        )
        .def("contains_species",
            &gridfire::reaction::ReactionSet::contains_species,
            py::arg("species"),
            "Check if any reaction in the set involves the given species."
        )
        .def(
            "contains_reactant",
            &gridfire::reaction::ReactionSet::contains_reactant,
            py::arg("species"),
            "Check if any reaction in the set has the species as a reactant."
        )
        .def(
            "contains_product",
            &gridfire::reaction::ReactionSet::contains_product,
            py::arg("species"),
            "Check if any reaction in the set has the species as a product."
        )
        .def(
            "__getitem__",
            py::overload_cast<size_t>(&gridfire::reaction::ReactionSet::operator[], py::const_),
            py::arg("index"),
            "Get a LogicalReaclibReaction by index."
        )
        .def(
            "__getitem___",
            py::overload_cast<const std::string_view&>(&gridfire::reaction::ReactionSet::operator[], py::const_),
            py::arg("id"),
            "Get a LogicalReaclibReaction by its ID."
        )
        .def(
            "__eq__",
            &gridfire::reaction::ReactionSet::operator==,
            py::arg("LogicalReactionSet"),
            "Equality operator for LogicalReactionSets based on their contents."
        )
        .def(
            "__ne__",
            &gridfire::reaction::ReactionSet::operator!=,
            py::arg("LogicalReactionSet"),
            "Inequality operator for LogicalReactionSets based on their contents."
        )
        .def(
            "hash",
            &gridfire::reaction::ReactionSet::hash,
            py::arg("seed") = 0,
            "Compute a hash for the LogicalReactionSet based on its contents."
        )
        .def(
            "__repr__",
            [](const gridfire::reaction::ReactionSet& self) {
                std::stringstream ss;
                ss << self;
                return ss.str();
            }
        )
        .def(
            "getReactionSetSpecies",
            &gridfire::reaction::ReactionSet::getReactionSetSpecies,
             "Get all species involved in the reactions of the set as a set of Species objects."
        );

        m.def(
            "packReactionSet",
            &gridfire::reaction::packReactionSet,
            py::arg("reactionSet"),
            "Convert a ReactionSet to a LogicalReactionSet by aggregating reactions with the same peName."
        );

        m.def(
            "get_all_reactions",
            &gridfire::reaclib::get_all_reaclib_reactions,
            "Get all reactions from the REACLIB database."
        );
}
