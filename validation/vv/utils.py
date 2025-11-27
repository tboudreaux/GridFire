from fourdst.composition import Composition
from fourdst.composition import CanonicalComposition
from fourdst.atomic import Species
from gridfire.type import NetIn

def rescale_composition(comp_ref : Composition, ZZs : float, Y_primordial : float = 0.248) -> Composition:
    CC : CanonicalComposition  = comp_ref.getCanonicalComposition()

    dY_dZ = (CC.Y - Y_primordial) / CC.Z

    Z_new = CC.Z * (10**ZZs)
    Y_bulk_new = Y_primordial + (dY_dZ * Z_new)
    X_new = 1.0 - Z_new - Y_bulk_new

    if X_new < 0: raise ValueError(f"ZZs={ZZs} yields unphysical composition (X < 0)")

    ratio_H  = X_new / CC.X if CC.X > 0 else 0
    ratio_He = Y_bulk_new / CC.Y if CC.Y > 0 else 0
    ratio_Z  = Z_new / CC.Z if CC.Z > 0 else 0

    Y_new_list = []
    newComp : Composition = Composition()
    s: Species
    for s  in comp_ref.getRegisteredSpecies():
        Xi_ref = comp_ref.getMassFraction(s)

        if s.el() == "H":
            Xi_new = Xi_ref * ratio_H
        elif s.el() == "He":
            Xi_new = Xi_ref * ratio_He
        else:
            Xi_new = Xi_ref * ratio_Z

        Y = Xi_new / s.mass()
        newComp.registerSpecies(s)
        newComp.setMolarAbundance(s, Y)

    return newComp

def init_composition(ZZs : float = 0) -> Composition:
    Y_solar = [7.0262E-01, 9.7479E-06, 6.8955E-02, 2.5000E-04, 7.8554E-05, 6.0144E-04, 8.1031E-05, 2.1513E-05]
    S = ["H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"]
    return rescale_composition(Composition(S, Y_solar), ZZs)


def init_netIn(temp: float, rho: float, time: float, comp: Composition) -> NetIn:
    n : NetIn = NetIn()
    n.temperature = temp
    n.density = rho
    n.tMax = time
    n.dt0 = 1e-12
    n.composition = comp
    return n

def years_to_seconds(years: float) -> float:
    return years * 3.1536e7
