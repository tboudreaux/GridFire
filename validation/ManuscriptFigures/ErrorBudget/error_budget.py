import numpy as np
import pandas as pd
from IPython.core.pylabtools import figsize
from gridfire.solver import PointSolver, PointSolverContext
from gridfire.policy import MainSequencePolicy

from gridfire.engine import GraphEngine, MultiscalePartitioningEngineView, AdaptiveEngineView
from gridfire.engine import NetworkBuildDepth

from scipy.signal import find_peaks

from gridfire.config import GridFireConfig

from fourdst.composition import Composition
from scipy.integrate import trapezoid

from fourdst.composition import CanonicalComposition
from fourdst.atomic import Species
from gridfire.type import NetIn, NetOut

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
    Y_solar = [7.0262E-01, 1.7479E-06, 6.8955E-02, 2.5000E-04, 7.8554E-05, 6.0144E-04, 8.1031E-05, 2.1513E-05]
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

def quantify_engine_error(df_base, df_approx, r_base: NetOut, r_approx: NetOut, species_list, floor_val=1e-30):
    temporal_results = {}
    final_state_results = {}

    t_base = df_base['t'].values

    tracking_cols = ['eps'] + species_list

    for col in tracking_cols:
        if col not in df_base.columns or col not in df_approx.columns:
            continue

        y_base = df_base[col].values

        interpolator = interp1d(
            df_approx['t'],
            df_approx[col],
            kind='linear',
            bounds_error=False,
            fill_value=(df_approx[col].iloc[0], df_approx[col].iloc[-1])
        )

        y_approx_interp = interpolator(t_base)

        abs_diff = np.abs(y_approx_interp - y_base)

        rel_diff = abs_diff / np.maximum(np.abs(y_base), floor_val)

        l2_diff = np.sqrt(trapezoid(abs_diff**2, x=t_base))
        l2_base = np.sqrt(trapezoid(y_base**2, x=t_base))

        temporal_results[col] = {
            'Max Rel Error (Temporal)': np.max(rel_diff),
            'L2 Rel Error (Temporal)': l2_diff / max(l2_base, floor_val)
        }
    def calc_rel_err(val_approx, val_base):
        return abs(val_approx - val_base) / max(abs(val_base), floor_val)

    final_state_results['Energy'] = {
        'Final Rel Error': calc_rel_err(r_approx.energy, r_base.energy)
    }

    final_state_results['Neutrino Loss'] = {
        'Final Rel Error': calc_rel_err(r_approx.specific_neutrino_energy_loss, r_base.specific_neutrino_energy_loss)
    }

    for sp in species_list:
        try:
            val_base = r_base.composition[sp]
            val_approx = r_approx.composition[sp]
            final_state_results[f"Final {sp}"] = {
                'Final Rel Error': calc_rel_err(val_approx, val_base)
            }
        except (KeyError, TypeError, AttributeError):
            pass

    return pd.DataFrame(temporal_results).T, pd.DataFrame(final_state_results).T

def main(save_show):
    C = init_composition()
    netIn = init_netIn(1.5e7, 1.6e2, years_to_seconds(10e9), C)


    stepLogger = StepLogger()

    engine_graph = GraphEngine(C, 4)

    solver_ctx_graph = PointSolverContext(engine_graph.constructStateBlob())
    solver_ctx_graph.stdout_logging = True
    solver_ctx_graph.callback = lambda ctx: stepLogger.log_step(ctx)

    solver_single = PointSolver(engine_graph)

    r_graph = solver_single.evaluate(solver_ctx_graph, netIn, False, False)
    df_graph : pd.DataFrame = stepLogger.df
    stepLogger.reset()

    QSE_engine = MultiscalePartitioningEngineView(engine_graph)
    solver_ctx_graph_qse = PointSolverContext(QSE_engine.constructStateBlob(engine_graph.constructStateBlob()))
    solver_ctx_graph_qse.stdout_logging = True
    solver_ctx_graph_qse.callback = lambda ctx: stepLogger.log_step(ctx)

    solver_QSE = PointSolver(QSE_engine)
    r_qse = solver_QSE.evaluate(solver_ctx_graph_qse, netIn, False, False)

    df_qse = stepLogger.df
    stepLogger.reset()

    # policy = MainSequencePolicy(C)
    # construct = policy.construct()
    # solver_AE_QSE = PointSolver(construct.engine)
    # solver_ctx_graph_qse_ae = PointSolverContext(construct.scratch_blob)
    # solver_ctx_graph_qse_ae.callback = lambda ctx: stepLogger.log_step(ctx)
    # solver_ctx_graph_qse_ae.stdout_logging = False
    #
    # r_ae_qse = solver_AE_QSE.evaluate(solver_ctx_graph_qse_ae, netIn, False, False)
    #
    # df_ae_qse = stepLogger.df
    # stepLogger.reset()

    # fig, axs = plt.subplots(2, 1, figsize=(10, 7))
    S = ["H-1", "He-4", "C-12", "N-14", "O-16", "Mg-24"]
    t = np.logspace(7, 17.5, 5000)
    # for spID, sp in enumerate(S):
    #     gf = interp1d(df_graph.t, df_graph[sp])
    #     qf = interp1d(df_qse.t, df_qse[sp])
    #
    #     ax = axs[0]
    #     ax.loglog(t, gf(t), 'o-', color=f"C{spID}")
    #     ax.loglog(t, qf(t), 'o', color=f"C{spID}", linestyle='dashed')
    #
    #     ax.text(1, df_graph[sp].iloc[0]*1.1, sp, fontsize=12, color=f"C{spID}")
    #
    #     ax = axs[1]
    #     ax.semilogx(t, (qf(t)-gf(t))/gf(t), color=f"C{spID}")
    #
    # axs[1].set_xlabel("Time [s]", fontsize=15)
    # axs[0].set_ylabel("Molar Abundance [mol/g]", fontsize=15)
    # axs[1].set_ylabel("Relative Error", fontsize=15)
    #
    # fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    # ge = interp1d(df_graph.t, df_graph.eps)
    # qe = interp1d(df_qse.t, df_qse.eps)
    # ax.loglog(t, np.abs((qe(t) - ge(t)) / ge(t)))

    temporal_err_qse, final_err_qse = quantify_engine_error(
        df_base=df_graph,
        df_approx=df_qse,
        r_base=r_graph,
        r_approx=r_qse,
        species_list=S
    )

    qse_rel_eps_error = (df_graph.eps.iloc[-1] - df_qse.eps.iloc[-1])/df_qse.eps.iloc[-1]


    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    # ax.semilogx(df_graph.t, df_graph["H-2"], 'o-', color='red')
    # ax.semilogx(df_qse.t, df_qse["H-2"], 'o', color='blue', linestyle='dashed')

    graph_h1 = interp1d(df_graph.t, df_graph["H-1"])
    qse_h1 = interp1d(df_qse.t, df_qse["H-1"])
    graph_h2 = interp1d(df_graph.t, df_graph["H-2"])
    qse_h2 = interp1d(df_qse.t, df_qse["H-2"])

    graph_DH = graph_h2(t)/graph_h1(t)
    qse_DH = qse_h2(t)/qse_h1(t)

    dex_diff = np.abs(np.log10(graph_h2(t)) - np.log10(qse_h2(t)))
    dex_dh_diff = np.abs(np.log10(graph_DH) - np.log10(qse_DH))
    # ax.semilogx(t, dex_diff, color='green')
    ax.loglog(t, dex_dh_diff, color='black')
    # ax.semilogx(t, qse_h2(t)/qse_h1(t), color='green')
    ax.set_xlabel("Time [s]", fontsize=17)
    ax.set_ylabel(r"$\left|\Delta\log_{10}\right|$ [dex]", fontsize=17)

    if save_show == ShowSave.SAVE:
        plt.savefig("DHErr.pdf")
        plt.close()
    else:
        plt.show()

    sums_qse = {}
    sums_graph = {}
    symbols = {}

    for sp, y in r_qse.composition:
        z = sp.z()
        symbols[z] = sp.el()

        y_graph = r_graph.composition.getMolarAbundance(sp)

        sums_qse[int(z)] = sums_qse.get(z, 0.0) + y
        sums_graph[int(z)] = sums_graph.get(z, 0.0) + y_graph

    print(sums_qse[3])
    print(sums_graph[3])

    z_list = sorted(sums_qse.keys())
    dex_list = []

    symbols = [val for key, val in symbols.items()]

    for z in z_list:
        total_qse = sums_qse[z]
        total_graph = sums_graph[z]

        if total_graph > 1e-13 and total_qse > 1e-13:
            offset = np.log10(total_qse / total_graph)
        else:
            if z >= 14:
                offset = np.nan 
            else:
                offset = 0.0

        dex_list.append(offset)
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    data = sorted(zip(z_list, symbols, dex_list), key=lambda x: x[0])
    sorted_z, sorted_symbols, sorted_dex = zip(*data)
    print(sorted_symbols)
    print(sorted_dex)
    # 2. Create the plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    print(sorted_symbols)
    bars = ax.bar(sorted_symbols, sorted_dex, color='grey', edgecolor='grey', alpha=0.8)

    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_xlabel('Element', fontsize=25)
    ax.set_ylabel('Offset [dex]', fontsize=25)

    if save_show == ShowSave.SAVE:
        plt.savefig("DexElementalOffset.pdf")
        plt.close()


    e_graph = interp1d(df_graph.t, df_graph.eps)
    e_qse = interp1d(df_qse.t, df_qse.eps)

    dex_eps_diff = np.log10(e_graph(t)) - np.log10(e_qse(t))
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    ax.semilogx(t, dex_eps_diff, color='black')
    ax.set_xlabel("Time [s]", fontsize=25)
    ax.set_ylabel("Offset [dex]", fontsize=25)

    if save_show == ShowSave.SAVE:
        plt.savefig("DexEpsOffset.pdf")
        plt.close()


    if save_show == ShowSave.SHOW:
        plt.show()

    print("=== QSE ===")
    print(temporal_err_qse)
    print(final_err_qse)
    print(f"Relative ε error: {qse_rel_eps_error}")

    print(f"Neutrino Loss Difference [dex]: {np.log10(r_graph.specific_neutrino_energy_loss) - np.log10(r_qse.specific_neutrino_energy_loss)}")

if __name__ == "__main__":
    import argparse
    app = argparse.ArgumentParser(prog="Derivative Smoothness", description="Generate of view plots of derivative smoothness")
    app.add_argument("-s", type=ShowSave, default=ShowSave.SHOW, choices=list(ShowSave), help="Whether to show or save the generated plot")

    args = app.parse_args()
    main(args.s)
