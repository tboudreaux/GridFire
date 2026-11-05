"""
GridFire type bindings
"""
from __future__ import annotations
import fourdst._phys.atomic
import fourdst._phys.composition
import typing
__all__: list[str] = ['InstantaneousRates', 'NetIn', 'NetOut', 'SharedCompositionResult', 'ShellState']
class InstantaneousRates:
    """
    Instantaneous rates at one state, with nothing integrated (see evaluateRates)
    """
    def __init__(self) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def dEps_dRho(self) -> float:
        """
        d(eps_nuc)/d(rho) at fixed T and composition (erg/g/s per g/cm^3)
        """
    @property
    def dEps_dT(self) -> float:
        """
        d(eps_nuc)/dT at fixed rho and composition (erg/g/s/K)
        """
    @property
    def dydt(self) -> dict[fourdst._phys.atomic.Species, float]:
        """
        dY/dt of every network species (mol/g/s)
        """
    @property
    def energy_channel_names(self) -> list[str]:
        """
        Names of the energy channels
        """
    @property
    def eps_by_channel(self) -> list[float]:
        """
        eps_nuc by energy channel (erg/g/s); sums to eps_nuc
        """
    @property
    def eps_nuc(self) -> float:
        """
        Specific nuclear energy generation rate (erg/g/s), net of reaction neutrino losses
        """
    @property
    def neutrino_energy_loss_rate(self) -> float:
        """
        Reaction neutrino energy loss rate (erg/g/s)
        """
    @property
    def neutrino_flux(self) -> float:
        """
        Reaction neutrino emission rate (neutrinos/g/s)
        """
    @property
    def reaction_energy_rate(self) -> list[float]:
        """
        Energy released per reaction, net of neutrinos (erg/g/s); sums to eps_nuc
        """
    @property
    def reaction_ids(self) -> list[str]:
        """
        Reaction identifiers, e.g. 'p(p,e+)d'
        """
    @property
    def reaction_molar_flow(self) -> list[float]:
        """
        Net molar flow of each reaction (mol/g/s)
        """
    @property
    def reaction_neutrino_energy_loss_rate(self) -> list[float]:
        """
        Neutrino energy loss per reaction (erg/g/s)
        """
    @property
    def reaction_neutrino_rate(self) -> list[float]:
        """
        Neutrinos emitted per reaction (1/g/s)
        """
class NetIn:
    composition: fourdst._phys.composition.Composition
    def __init__(self) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def density(self) -> float:
        ...
    @density.setter
    def density(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def dt0(self) -> float:
        ...
    @dt0.setter
    def dt0(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def energy(self) -> float:
        ...
    @energy.setter
    def energy(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def tMax(self) -> float:
        ...
    @tMax.setter
    def tMax(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def temperature(self) -> float:
        ...
    @temperature.setter
    def temperature(self, arg0: typing.SupportsFloat) -> None:
        ...
class NetOut:
    def __repr__(self) -> str:
        ...
    @property
    def converged(self) -> bool:
        """
        False if this zone could not be integrated to tMax (multi-zone solvers); its other values are then zero
        """
    @property
    def failure_message(self) -> str:
        """
        Why the zone failed (empty if converged)
        """
    @property
    def time_reached(self) -> float:
        """
        Time integrated to (s): tMax on success, the last accepted step on failure
        """
    @property
    def composition(self) -> fourdst._phys.composition.Composition:
        ...
    @property
    def dEps_dRho(self) -> float:
        ...
    @property
    def dEps_dT(self) -> float:
        ...
    @property
    def dlnEps_dlnRho(self) -> float:
        """
        Logarithmic derivative dln(eps_nuc)/dln(rho) at the final state (0 if eps_nuc == 0)
        """
    @property
    def dlnEps_dlnT(self) -> float:
        """
        Logarithmic derivative dln(eps_nuc)/dln(T) at the final state (0 if eps_nuc == 0)
        """
    @property
    def energy(self) -> float:
        ...
    @property
    def energy_channel_names(self) -> list[str]:
        """
        Names of the energy channels (e.g. H_pp, H_CNO, He, ...)
        """
    @property
    def eps_by_channel(self) -> list[float]:
        """
        Instantaneous eps_nuc (erg/g/s) attributed to each energy channel; sums to eps_nuc
        """
    @property
    def eps_channels(self) -> dict[str, float]:
        """
        Mapping from energy channel name to its instantaneous eps_nuc contribution (erg/g/s)
        """
    @property
    def eps_nuc(self) -> float:
        """
        Instantaneous specific nuclear energy generation rate (erg/g/s) at the final state, net of reaction neutrino losses
        """
    @property
    def neutrino_energy_loss_rate(self) -> float:
        """
        Instantaneous reaction neutrino energy loss rate (erg/g/s) at the final state
        """
    @property
    def num_steps(self) -> int:
        ...
    @property
    def specific_neutrino_energy_loss(self) -> float:
        ...
    @property
    def specific_neutrino_flux(self) -> float:
        ...
class SharedCompositionResult:
    """
    Result of burning one composition shared by several shells
    """
    @property
    def group(self) -> NetOut:
        """
        The shared composition, step count and time reached; energy is the mass-weighted mean
        """
    @property
    def shell_energy(self) -> list[float]:
        """
        Specific energy released in each shell (erg/g), net of neutrinos
        """
    @property
    def shell_neutrino_energy_loss(self) -> list[float]:
        """
        Reaction neutrino energy lost from each shell (erg/g)
        """
    @property
    def shell_neutrino_flux(self) -> list[float]:
        """
        Reaction neutrinos emitted per gram of each shell
        """
class ShellState:
    """
    One shell of a shared-composition group: its state and its share of the group's mass
    """
    density: float
    temperature: float
    weight: float
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, temperature: typing.SupportsFloat, density: typing.SupportsFloat, weight: typing.SupportsFloat) -> None:
        ...
