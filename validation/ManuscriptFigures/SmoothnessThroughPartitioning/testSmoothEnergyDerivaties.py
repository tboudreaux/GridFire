import numpy as np
from IPython.core.pylabtools import figsize
from gridfire.solver import PointSolver, PointSolverContext
from gridfire.policy import MainSequencePolicy

from scipy.signal import find_peaks

from gridfire.config import GridFireConfig

from fourdst.composition import Composition
from scipy.integrate import trapezoid

from fourdst.composition import CanonicalComposition
from fourdst.atomic import Species
from gridfire.type import NetIn

import matplotlib.pyplot as plt


## Note that my default style uses tex rendering. If you do not have tex installed
## simply comment out this line
plt.style.use("../utils/pub.mplstyle")

from scipy.interpolate import interp1d, CubicSpline

from enum import Enum

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../utils")))

from logger import StepLogger

class ShowSave(Enum):
    SHOW="SHOW"
    SAVE="SAVE"

    def __str__(self):
        return self.value

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

def main(save_show):
    C = init_composition()
    netIn = init_netIn(1.5e7, 160, years_to_seconds(10e9), C)
    policy = MainSequencePolicy(C)
    construct = policy.construct()

    # 3e-8 and 1e-24 are the default tolerances we adopt as testing indicates it works well for
    # main sequence evolution. We encorage researchers to trial various relative and
    # absolute thresholds
    # config = GridFireConfig()
    # config.solver.pointSolver.trigger.boundaryFlux.relativeThreshold = 3e-8
    # config.solver.pointSolver.trigger.boundaryFlux.absoluteThreshold = 1e-24
    # solver = PointSolver(construct.engine, config)

    solver = PointSolver(construct.engine)
    solver_ctx = PointSolverContext(construct.scratch_blob)

    stepLogger = StepLogger()
    solver_ctx.callback = lambda ctx: stepLogger.log_step(ctx);
    solver.evaluate(solver_ctx, netIn, False, False)

    df = stepLogger.df


    fig, axs = plt.subplots(2, 1, figsize=(17, 10))

    t = np.linspace(df.t.min(), df.t.max(), 1000)

    # Note we are not plotting Ne-20 as its molar abundance is so close to N-14 that it makes it hard to
    # distinguish that species
    PlottingSpecies = ["H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Mg-24"]
    stable_index = 10

    for sp in PlottingSpecies:
        x = df.t[stable_index:]
        y = df[sp][stable_index:]
        axs[0].loglog(x, y)
        axs[1].semilogx(x, np.gradient(y, x))

        axs[0].text(x.iloc[0], y.iloc[0]*1.1, sp, fontsize=12)

    axs[0].set_ylabel("$Y$ [mol/g]", fontsize=23)
    axs[1].set_ylabel(r"$\frac{dY}{dt}$ [mol/g/s]", fontsize=23)
    axs[1].set_xlabel("Time [s]")

    ax_eps = axs[0].twinx()
    ax_deps = axs[1].twinx()

    ax_eps.set_ylabel(r"$\epsilon$ [erg/g/s]", rotation=270, labelpad=25, fontsize=23)
    ax_deps.set_ylabel(r"$\frac{d\epsilon}{dt}$ [erg/g/s$^2$]", rotation=270, labelpad=25, fontsize=23)

    ax_eps.axvline(1.008e+15, color='grey', linestyle='dashed')
    ax_deps.axvline(1.008e+15, color='grey', linestyle='dashed')

    ax_eps.loglog(df.t[stable_index:], df.eps[stable_index:], color='red', linestyle='dashed')
    ax_eps.text(df.t[stable_index:].iloc[0]*1.05, df.eps[stable_index:].iloc[0]*3, r"$\epsilon$", rotation=25, fontsize=20)

    ax_deps.semilogx(df.t[stable_index:], np.gradient(df.eps[stable_index:], df.t[stable_index:]), color='red', linestyle='dashed')


    if save_show == ShowSave.SHOW:
        plt.show()
    else:
        plt.savefig("smoothness_plot.pdf")
    plt.close()


    t = df.t.values
    eps = df.eps.values

    # Use this plot to determine the index to test removal of
    # fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    # ax.plot(np.gradient(eps, t))
    # ax.grid()
    # plt.show()

    idx = 156
    t1 = t
    eps1 = eps

    t2 = np.delete(t, idx)
    eps2 = np.delete(eps, idx)

    f_deps_1 = interp1d(t1, np.gradient(eps1, t1))
    f_deps_2 = interp1d(t2, np.gradient(eps2, t2))

    int_deps_1 = trapezoid(f_deps_1(t), t)
    int_deps_2 = trapezoid(f_deps_2(t), t)

    rel_err = (int_deps_1 - int_deps_2) / int_deps_2
    print(f"Rel Error: {rel_err:+0.3E}")

    window = 10
    indices = np.arange(idx - window, idx + window + 1)
    indices_no_gap = np.delete(indices, window)

    clean_t = t[indices_no_gap]
    clean_eps = eps[indices_no_gap]
    spline = CubicSpline(clean_t, clean_eps)

    eps_predicted = spline(t[idx])
    eps_actual = eps[idx]

    absolute_jump = np.abs(eps_actual - eps_predicted)
    relative_jump = absolute_jump / eps_actual

    print(f"Local Discontinuity at index {idx}: {relative_jump:.3%}")

    E_actual = trapezoid(eps, t)

    t_clean = np.delete(t, idx)
    eps_clean_points = np.delete(eps, idx)
    spline = CubicSpline(t_clean, eps_clean_points)

    eps_smooth = np.copy(eps)
    eps_smooth[idx] = spline(t[idx])

    E_smooth = trapezoid(eps_smooth, t)

    total_rel_error = (E_actual - E_smooth) / E_smooth
    print(f"Total Relative Energy Error: {total_rel_error:+0.12E}")







if __name__ == "__main__":
    import argparse
    app = argparse.ArgumentParser(prog="Derivative Smoothness", description="Generate of view plots of derivative smoothness")
    app.add_argument("-s", type=ShowSave, default=ShowSave.SHOW, choices=list(ShowSave), help="Whether to show or save the generated plot")

    args = app.parse_args()
    main(args.s)
