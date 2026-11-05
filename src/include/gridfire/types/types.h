/* ***********************************************************************
//
//   Copyright (C) 2025 -- The 4D-STAR Collaboration
//   File Authors: Emily Boudreaux, Aaron Dotter
//   Last Modified: Nov 24, 2025
//
//   GridFire is free software; you can use it and/or modify
//   it under the terms and restrictions the GNU General Library Public
//   License version 3 (GPLv3) as published by the Free Software Foundation.
//
//   GridFire is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Library General Public License for more details.
//
//   You should have received a copy of the GNU Library General Public License
//   along with this software; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// *********************************************************************** */
#pragma once

#include "fourdst/composition/composition.h"
#include "cppad/cppad.hpp"
#include <format>
#include <string>
#include <map>
#include <vector>

namespace gridfire {
    struct NetIn {
        fourdst::composition::Composition composition; ///< Composition of the network
        double tMax; ///< Maximum time
        double dt0; ///< Initial time step
        double temperature; ///< Temperature in Kelvin
        double density; ///< Density in g/cm^3
        double energy; ///< Energy in ergs
        double culling = 0.0; ///< Culling threshold for reactions (default is 0.0, meaning no culling)
    };

    struct NetOut {
        fourdst::composition::Composition composition; ///< Composition of the network after evaluation (the input composition for a zone that failed)
        int num_steps = 0; ///< Number of steps taken in the evaluation
        double energy = 0.0; ///< Energy in ergs after evaluation
        double dEps_dT = 0.0; ///< Partial derivative of energy generation rate with respect to temperature
        double dEps_dRho = 0.0; ///< Partial derivative of energy generation rate with respect to density
        double specific_neutrino_energy_loss = 0.0; ///< Neutrino energy loss rate in ergs/g/s
        double specific_neutrino_flux = 0.0; ///< Total neutrino flux in neutrinos/g/s

        // --- Outcome (multi-zone solvers report a failed zone here; single-zone solvers throw instead) ---
        bool converged = true; ///< false if this zone could not be integrated to tMax; its other values are then zero
        double time_reached = 0.0; ///< Time integrated to [s]: tMax on success, the last accepted step on failure
        std::string failure_message; ///< Why the zone failed (empty if converged)

        // --- Instantaneous quantities at the final state (same composition, T, and rho as dEps_dT and dEps_dRho) ---
        double eps_nuc = 0.0; ///< Instantaneous specific nuclear energy generation rate, net of reaction neutrino losses, in erg/g/s
        double dlnEps_dlnT = 0.0; ///< Logarithmic derivative (T/eps_nuc) * dEps_dT at fixed rho and composition (0 if eps_nuc == 0)
        double dlnEps_dlnRho = 0.0; ///< Logarithmic derivative (rho/eps_nuc) * dEps_dRho at fixed T and composition (0 if eps_nuc == 0)
        double neutrino_energy_loss_rate = 0.0; ///< Instantaneous reaction neutrino energy loss rate in erg/g/s
        std::vector<std::string> energy_channel_names; ///< Names of the energy channels (e.g. H_pp, H_CNO, He, ...)
        std::vector<double> eps_by_channel; ///< Instantaneous eps_nuc attributed to each channel in erg/g/s; sums to eps_nuc

        friend std::ostream& operator<<(std::ostream& os, const NetOut& netOut) {
            os << "NetOut(composition=" << netOut.composition << ", num_steps=" << netOut.num_steps << ", ε=" << netOut.energy << ", dε/dT=" << netOut.dEps_dT << ", dε/dρ=" << netOut.dEps_dRho << ", ε_nuc=" << netOut.eps_nuc << ")";
            return os;
        }
    };

    /**
     * @brief Instantaneous rates at one state (composition, T, rho), with nothing integrated.
     *
     * Returned by SingleZoneNetworkSolver::evaluateRates and MultiZoneNetworkSolver::evaluateRates, and what the C
     * API's gf_evaluate reports. eps_nuc, dEps_dT and dEps_dRho match the final-state values evolve() reports for the
     * same state.
     */
    struct InstantaneousRates {
        double eps_nuc = 0.0; ///< Specific nuclear energy generation rate, net of reaction neutrino losses, in erg/g/s
        double dEps_dT = 0.0; ///< d(eps_nuc)/dT at fixed rho and composition, in erg/g/s/K
        double dEps_dRho = 0.0; ///< d(eps_nuc)/d(rho) at fixed T and composition, in erg/g/s/(g/cm^3)
        double neutrino_energy_loss_rate = 0.0; ///< Reaction neutrino energy loss rate in erg/g/s
        double neutrino_flux = 0.0; ///< Reaction neutrino emission rate in neutrinos/g/s
        std::map<fourdst::atomic::Species, double> dydt; ///< dY/dt of every network species, in mol/g/s
        std::vector<std::string> energy_channel_names; ///< Names of the energy channels (empty if not requested)
        std::vector<double> eps_by_channel; ///< eps_nuc by energy channel, in erg/g/s; sums to eps_nuc (empty if not requested)
    };

    /**
     * @brief Concept for types allowed in engine calculations.
     *
     * This concept restricts template parameters to either double or CppAD::AD<double>,
     * enabling both standard and automatic differentiation types.
     */
    template<typename T>
    concept IsArithmeticOrAD = std::is_same_v<T, double> || std::is_same_v<T, CppAD::AD<double>>;

} // namespace nuclearNetwork

template<>
struct std::formatter<gridfire::NetIn> : std::formatter<std::string> {
    auto format(const gridfire::NetIn& netIn, auto& ctx) {
        std::string output = "NetIn(, tMax=" + std::to_string(netIn.tMax) +
                             ", dt0=" + std::to_string(netIn.dt0) +
                             ", temperature=" + std::to_string(netIn.temperature) +
                             ", density=" + std::to_string(netIn.density) +
                             ", energy=" + std::to_string(netIn.energy) + ")";
        return std::formatter<std::string>::format(output, ctx);
    }
};

template <>
struct std::formatter<gridfire::NetOut> : std::formatter<std::string> {
    auto format(const gridfire::NetOut& netOut, auto& ctx) {
        std::string output = "NetOut(, num_steps=" + std::to_string(netOut.num_steps) +
                             ", energy=" + std::to_string(netOut.energy) +
                             ", dEps_dT=" + std::to_string(netOut.dEps_dT) +
                             ", dEps_dRho=" + std::to_string(netOut.dEps_dRho) +
                             ", eps_nuc=" + std::to_string(netOut.eps_nuc) + ")";
        return std::formatter<std::string>::format(output, ctx);
    }
};
