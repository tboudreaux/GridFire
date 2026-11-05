#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly
#include <format>

#include "bindings.h"

namespace py = pybind11;

#include "gridfire/types/types.h"

void register_type_bindings(const pybind11::module &m) {
    py::class_<gridfire::NetIn>(m, "NetIn")
    .def(py::init<>())
    .def_readwrite("composition", &gridfire::NetIn::composition)
    .def_readwrite("tMax", &gridfire::NetIn::tMax)
    .def_readwrite("dt0", &gridfire::NetIn::dt0)
    .def_readwrite("temperature", &gridfire::NetIn::temperature)
    .def_readwrite("density", &gridfire::NetIn::density)
    .def_readwrite("energy", &gridfire::NetIn::energy)
    .def("__repr__", [](const gridfire::NetIn &netIn) {
        std::stringstream ss;
        ss << "NetIn(composition=" << netIn.composition
           << ", tMax=" << netIn.tMax
           << ", dt0=" << netIn.dt0
           << ", temperature=" << netIn.temperature
           << ", density=" << netIn.density
           << ", energy=" << netIn.energy << ")";
        return ss.str();
    });

    py::class_<gridfire::ShellState>(m, "ShellState", "One shell of a shared-composition group: its state and its share of the group's mass")
        .def(py::init<>())
        .def(py::init([](const double temperature, const double density, const double weight) {
            return gridfire::ShellState{temperature, density, weight};
        }), py::arg("temperature"), py::arg("density"), py::arg("weight"))
        .def_readwrite("temperature", &gridfire::ShellState::temperature, "Temperature (K)")
        .def_readwrite("density", &gridfire::ShellState::density, "Density (g/cm^3)")
        .def_readwrite("weight", &gridfire::ShellState::weight, "Mass fraction of the group in this shell");

    py::class_<gridfire::SharedCompositionResult>(m, "SharedCompositionResult", "Result of burning one composition shared by several shells")
        .def_readonly("group", &gridfire::SharedCompositionResult::group, "The shared composition, step count and time reached; energy is the mass-weighted mean")
        .def_readonly("shell_energy", &gridfire::SharedCompositionResult::shell_energy, "Specific energy released in each shell (erg/g), net of neutrinos")
        .def_readonly("shell_neutrino_energy_loss", &gridfire::SharedCompositionResult::shell_neutrino_energy_loss, "Reaction neutrino energy lost from each shell (erg/g)")
        .def_readonly("shell_neutrino_flux", &gridfire::SharedCompositionResult::shell_neutrino_flux, "Reaction neutrinos emitted per gram of each shell");

    py::class_<gridfire::InstantaneousRates>(m, "InstantaneousRates", "Instantaneous rates at one state, with nothing integrated (see evaluateRates)")
        .def(py::init<>())
        .def_readonly("eps_nuc", &gridfire::InstantaneousRates::eps_nuc, "Specific nuclear energy generation rate (erg/g/s), net of reaction neutrino losses")
        .def_readonly("dEps_dT", &gridfire::InstantaneousRates::dEps_dT, "d(eps_nuc)/dT at fixed rho and composition (erg/g/s/K)")
        .def_readonly("dEps_dRho", &gridfire::InstantaneousRates::dEps_dRho, "d(eps_nuc)/d(rho) at fixed T and composition (erg/g/s per g/cm^3)")
        .def_readonly("neutrino_energy_loss_rate", &gridfire::InstantaneousRates::neutrino_energy_loss_rate, "Reaction neutrino energy loss rate (erg/g/s)")
        .def_readonly("neutrino_flux", &gridfire::InstantaneousRates::neutrino_flux, "Reaction neutrino emission rate (neutrinos/g/s)")
        .def_readonly("dydt", &gridfire::InstantaneousRates::dydt, "dY/dt of every network species (mol/g/s)")
        .def_readonly("energy_channel_names", &gridfire::InstantaneousRates::energy_channel_names, "Names of the energy channels")
        .def_readonly("eps_by_channel", &gridfire::InstantaneousRates::eps_by_channel, "eps_nuc by energy channel (erg/g/s); sums to eps_nuc")
        .def("__repr__", [](const gridfire::InstantaneousRates& rates) {
            return std::format("InstantaneousRates(ε_nuc = {}, dε/dT = {}, dε/dρ = {})", rates.eps_nuc, rates.dEps_dT, rates.dEps_dRho);
        });

    py::class_<gridfire::NetOut>(m, "NetOut")
        .def_readonly("composition", &gridfire::NetOut::composition)
        .def_readonly("num_steps", &gridfire::NetOut::num_steps)
        .def_readonly("energy", &gridfire::NetOut::energy)
        .def_readonly("dEps_dT", &gridfire::NetOut::dEps_dT)
        .def_readonly("dEps_dRho", &gridfire::NetOut::dEps_dRho)
        .def_readonly("specific_neutrino_flux", &gridfire::NetOut::specific_neutrino_flux)
        .def_readonly("specific_neutrino_energy_loss", &gridfire::NetOut::specific_neutrino_energy_loss)
        .def_readonly("converged", &gridfire::NetOut::converged, "False if this zone could not be integrated to tMax (multi-zone solvers); its other values are then zero")
        .def_readonly("time_reached", &gridfire::NetOut::time_reached, "Time integrated to (s): tMax on success, the last accepted step on failure")
        .def_readonly("failure_message", &gridfire::NetOut::failure_message, "Why the zone failed (empty if converged)")
        .def_readonly("eps_nuc", &gridfire::NetOut::eps_nuc, "Instantaneous specific nuclear energy generation rate (erg/g/s) at the final state, net of reaction neutrino losses")
        .def_readonly("dlnEps_dlnT", &gridfire::NetOut::dlnEps_dlnT, "Logarithmic derivative dln(eps_nuc)/dln(T) at the final state (0 if eps_nuc == 0)")
        .def_readonly("dlnEps_dlnRho", &gridfire::NetOut::dlnEps_dlnRho, "Logarithmic derivative dln(eps_nuc)/dln(rho) at the final state (0 if eps_nuc == 0)")
        .def_readonly("neutrino_energy_loss_rate", &gridfire::NetOut::neutrino_energy_loss_rate, "Instantaneous reaction neutrino energy loss rate (erg/g/s) at the final state")
        .def_readonly("energy_channel_names", &gridfire::NetOut::energy_channel_names, "Names of the energy channels (e.g. H_pp, H_CNO, He, ...)")
        .def_readonly("eps_by_channel", &gridfire::NetOut::eps_by_channel, "Instantaneous eps_nuc (erg/g/s) attributed to each energy channel; sums to eps_nuc")
        .def_property_readonly("eps_channels", [](const gridfire::NetOut& netOut) {
            py::dict channels;
            for (size_t i = 0; i < netOut.energy_channel_names.size() && i < netOut.eps_by_channel.size(); ++i) {
                channels[py::str(netOut.energy_channel_names[i])] = netOut.eps_by_channel[i];
            }
            return channels;
        }, "Mapping from energy channel name to its instantaneous eps_nuc contribution (erg/g/s)")
        .def("__repr__", [](const gridfire::NetOut &netOut) {
            std::string repr = std::format(
                "NetOut(<μ> = {} steps = {}, ε = {}, dε/dT = {}, dε/dρ = {}, ε_nuc = {})",
                netOut.composition.getMeanParticleMass(),
                netOut.num_steps,
                netOut.energy,
                netOut.dEps_dT,
                netOut.dEps_dRho,
                netOut.eps_nuc
            );
            return repr;
        });

}
