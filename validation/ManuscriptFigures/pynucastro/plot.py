import argparse
import json
import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from astropy import constants as const
from astropy import units as u
from enum import Enum
from typing import List, Dict, Any, Tuple
from fourdst.atomic import species as spdict

EXTERNAL_STYLE_PATH = "../ManuscriptFigures/utils/pub.mplstyle"

class PlotVariable(Enum):
    COMPOSITION = "composition"
    EPS = "eps"
    DT = "dt"

class OutputFormat(Enum):
    INTERACTIVE = "interactive"
    PDF = "pdf"
    PNG = "png"
    JPEG = "jpeg"

def discover_runs(base_dir: str) -> List[str]:
    runs = set()
    gf_ok_dir = os.path.join(base_dir, "GridFire", "Ok")

    if os.path.exists(gf_ok_dir):
        for fname in os.listdir(gf_ok_dir):
            if fname.endswith("_OKAY.json"):
                runs.add(fname.replace("_OKAY.json", ""))

    return sorted(list(runs))

def load_run_data(base_dir: str, run_name: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    gf_data = {}
    pynuc_data = {}

    gf_ok_path = os.path.join(base_dir, "GridFire", "Ok", f"{run_name}_OKAY.json")
    gf_err_path = os.path.join(base_dir, "GridFire", "Err", f"{run_name}_FAIL.json")
    pynuc_path = os.path.join(base_dir, "pynucastro", f"{run_name}_pynucastro.json")

    if os.path.exists(gf_ok_path):
        with open(gf_ok_path, 'r') as f:
            gf_data = json.load(f)
    elif os.path.exists(gf_err_path):
        with open(gf_err_path, 'r') as f:
            gf_data = json.load(f)

    if os.path.exists(pynuc_path):
        with open(pynuc_path, 'r') as f:
            pynuc_data = json.load(f)

    return gf_data, pynuc_data


def get_pynuc_eps(steps: List[Dict[str, Any]]) -> Tuple[np.ndarray, np.ndarray]:
    c_sq = (const.c.cgs.value)**2
    Na = const.N_A.cgs.value
    amu_to_g = const.u.cgs.value

    eps_history = []
    time_history = []
    last_Y = {}
    last_t = None

    for step in steps:
        t = step["t"]
        current_Y = step["Composition"]

        if last_t is None:
            last_t = t
            last_Y = current_Y.copy()
            eps_history.append(0.0)
            time_history.append(t)
            continue

        dt = t - last_t

        if dt > 0:
            dm_dt = 0.0

            all_species = set(current_Y.keys()) | set(last_Y.keys())

            for sp in all_species:
                y_curr = current_Y.get(sp, 0.0)
                y_prev = last_Y.get(sp, 0.0)
                dy = y_curr - y_prev
                if sp in spdict:
                    mass_g = spdict[sp].mass() * amu_to_g

                dm_dt += mass_g * (dy / dt)
            rate = -dm_dt * Na * c_sq
            eps_history.append(rate)
            time_history.append(t)

        last_t = t
        last_Y = current_Y.copy()

    return np.array(time_history), np.array(eps_history)

def _setup_axes(ax_main: plt.Axes, ax_diff: plt.Axes, var: PlotVariable, fig_opts: dict):
    ax_diff.set_xlabel("Time (s)")
    ax_main.set_xscale(fig_opts['x_scale'])
    ax_diff.set_xscale(fig_opts['x_scale'])

    ax_main.set_yscale(fig_opts['y_scale'])

    if var == PlotVariable.EPS:
        ax_main.set_ylabel("Cumulative Energy (eps)")
    elif var == PlotVariable.DT:
        ax_main.set_ylabel("Timestep Size (dt)")
    elif var == PlotVariable.COMPOSITION:
        ax_main.set_ylabel("Mass Fraction (X_i)")

    ax_diff.set_ylabel(r"$\Delta \log_{10}$")

def _plot_single_run(ax_main: plt.Axes, ax_diff: plt.Axes, run_name: str, var: PlotVariable,
                     base_dir: str, compare_pynuc: bool):
    gf_data, pynuc_data = load_run_data(base_dir, run_name)

    if not gf_data or gf_data.get("Metadata", {}).get("Status") == "Error":
        return

    gf_steps = gf_data.get("Steps", [])
    if not gf_steps:
        return

    if compare_pynuc and not pynuc_data:
        print(f"Warning: PyNucastro comparison requested but data not found for '{run_name}'.")

    t_gf = np.array([s["t"] for s in gf_steps])

    if var == PlotVariable.COMPOSITION:
        final_comp = gf_steps[-1]["Composition"]
        top_species = sorted(final_comp, key=final_comp.get, reverse=True)[:3]

        for spec in top_species:
            y_gf = np.array([s["Composition"].get(spec, 1e-30) for s in gf_steps])
            line, = ax_main.plot(t_gf, y_gf, label=f"{run_name} {spec} (GF)")

            if compare_pynuc and pynuc_data:
                pynuc_steps = pynuc_data.get("Steps", [])
                if not pynuc_steps:
                    continue

                t_pynuc = np.array([s["t"] for s in pynuc_steps])
                y_pynuc = np.array([s["Composition"].get(spec, 1e-30) for s in pynuc_steps])

                ax_main.plot(t_pynuc, y_pynuc, '--', color=line.get_color(), label=f"{run_name} {spec} (PyNuc)")

                if len(t_pynuc) > 1:
                    f_interp = interp1d(t_pynuc, y_pynuc, kind='linear', bounds_error=False, fill_value=(y_pynuc[0], y_pynuc[-1]))
                    y_pynuc_interp = f_interp(t_gf)

                    log_diff = np.abs(np.log10(np.maximum(y_gf, 1e-30)) - np.log10(np.maximum(y_pynuc_interp, 1e-30)))
                    ax_diff.plot(t_gf, log_diff, color=line.get_color(), linestyle=':', label=f"Δ {spec}")

    elif var == PlotVariable.EPS:
        y_gf = np.array([s["eps"] for s in gf_steps])
        line, = ax_main.plot(t_gf, y_gf, label=f"{run_name} (GF)")

        if compare_pynuc and pynuc_data:
            pynuc_steps = pynuc_data.get("Steps", [])
            if pynuc_steps:
                s_t, s_e = get_pynuc_eps(pynuc_steps)

                if len(s_t) > 1:
                    s_cumE = cumulative_trapezoid(s_e, s_t, initial=0)

                    ax_main.plot(s_t, s_cumE, '--', color=line.get_color(), label=f"{run_name} (PyNuc)")

                    f_pynuc_interp = interp1d(s_t[np.isfinite(s_cumE)], s_cumE[np.isfinite(s_cumE)])
                    f_gf_interp = interp1d(t_gf, y_gf)


                    t_safe = np.logspace(
                        8,
                        np.log10(min(s_t.max(), t_gf.max())),
                        1000
                    )
                    y_pynuc_interp = f_pynuc_interp(t_safe)
                    y_gf_interp = f_gf_interp(t_safe)

                    pynuc_safe = np.maximum(np.abs(y_pynuc_interp), 1e-30)
                    gf_safe = np.maximum(np.abs(y_gf_interp), 1e-30)
                    log_diff = np.log10(gf_safe) - np.log10(pynuc_safe)
                    ax_diff.plot(t_safe, log_diff, color=line.get_color(), linestyle=':', label=f"Δ eps")
                    ax_main.set_xlim(1e8)

    elif var == PlotVariable.DT:
        y_gf = np.array([s["dt"] for s in gf_steps])
        ax_main.plot(t_gf, y_gf, label=f"{run_name} (GF)")

def _finalize_plot(fig: plt.Figure, ax_main: plt.Axes, ax_diff: plt.Axes, format_opt: OutputFormat, filename_base: str, is_subfigure: bool = False):
    ax_main.legend(loc='best', fontsize='small')

    if len(ax_diff.lines) > 0:
        ax_diff.legend(loc='best', fontsize='x-small')

    if is_subfigure:
        return

    fig.tight_layout()

    if format_opt != OutputFormat.INTERACTIVE:
        out_name = f"{filename_base}.{format_opt.value}"
        fig.savefig(out_name, format=format_opt.value, bbox_inches='tight')
        print(f"Saved figure: {out_name}")
        plt.close(fig)

def plot_data(runs: List[str], plot_vars: List[PlotVariable], base_dir: str,
              compare_pynuc: bool, format_opt: OutputFormat, fig_opts: dict):
    if not runs:
        print("No valid runs to plot.")
        return

    if fig_opts['use_ext_style']:
        try:
            plt.style.use(EXTERNAL_STYLE_PATH)
            print(f"Using external style sheet: {EXTERNAL_STYLE_PATH}")
        except Exception as e:
            print(f"Warning: Failed to load external style sheet. Error: {e}")
    elif fig_opts['style']:
        try:
            plt.style.use(fig_opts['style'])
        except OSError:
            print(f"Warning: Style '{fig_opts['style']}' not found. Using default.")

    plt.rcParams["figure.figsize"] = fig_opts['figsize']
    plt.rcParams["figure.dpi"] = fig_opts['dpi']

    for var in plot_vars:
        if fig_opts['merge_runs']:
            fig, (ax_main, ax_diff) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
            ax_main.set_title(f"Comparison of {var.value.upper()} (Merged Runs)")
            _setup_axes(ax_main, ax_diff, var, fig_opts)

            for run_name in runs:
                _plot_single_run(ax_main, ax_diff, run_name, var, base_dir, compare_pynuc)

            _finalize_plot(fig, ax_main, ax_diff, format_opt, f"ValidationPlot_Merged_{var.value}")

        else:
            num_runs = len(runs)
            cols = math.ceil(math.sqrt(num_runs))
            rows = math.ceil(num_runs / cols)

            base_w, base_h = fig_opts['figsize']
            fig = plt.figure(figsize=(base_w * cols, base_h * rows), layout='constrained')
            fig.suptitle(f"{var.value.upper()} Comparison", fontsize=16, fontweight='bold')

            subfigs_raw = fig.subfigures(rows, cols)
            if hasattr(subfigs_raw, 'flatten'):
                subfigs = subfigs_raw.flatten()
            else:
                subfigs = [subfigs_raw]

            for i, run_name in enumerate(runs):
                subfig = subfigs[i]
                subfig.suptitle(f"{run_name}", fontsize=12)

                axes = subfig.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
                ax_main, ax_diff = axes[0], axes[1]

                _setup_axes(ax_main, ax_diff, var, fig_opts)
                _plot_single_run(ax_main, ax_diff, run_name, var, base_dir, compare_pynuc)

                _finalize_plot(fig, ax_main, ax_diff, format_opt, "", is_subfigure=True)

            for j in range(num_runs, len(subfigs)):
                subfigs[j].set_visible(False)

            if format_opt != OutputFormat.INTERACTIVE:
                out_name = f"ValidationPlot_Grid_{var.value}.{format_opt.value}"
                fig.savefig(out_name, format=format_opt.value, bbox_inches='tight')
                print(f"Saved grid figure: {out_name}")
                plt.close(fig)

    if format_opt == OutputFormat.INTERACTIVE:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="GridFire Validation Suite Output Parser and Plotter")

    parser.add_argument("-d", "--data-dir", type=str, default="GF_Validation_Output",
                        help="Path to the directory containing the JSON output folders.")
    parser.add_argument("--runs", nargs="+", type=str, required=False,
                        help="Which validation runs to analyze. Use 'all' to process all available runs.")
    parser.add_argument("--plot", nargs="*", type=lambda x: PlotVariable[x.upper()], choices=list(PlotVariable), default=[],
                        help="Variables to plot. Leave empty to skip plotting.")
    parser.add_argument("--compare-pynucastro", action="store_true",
                        help="Include pynucastro data and calculate log residuals.")
    parser.add_argument("--merge-runs", action="store_true",
                        help="Merge all specified runs onto a single figure per variable. (Default: Grid layout of subfigures)")
    parser.add_argument("--x-scale", type=str, choices=["log", "linear"], default="log",
                        help="Scale for the x-axis (time). Default is 'log'.")
    parser.add_argument("--y-scale", type=str, choices=["log", "linear"], default="log",
                        help="Scale for the y-axis (main plots). Default is 'log'.")
    parser.add_argument("--format", type=lambda x: OutputFormat[x.upper()], choices=list(OutputFormat), default=OutputFormat.INTERACTIVE,
                        help="Output format for the plots. Default is interactive window.")
    parser.add_argument("--use-external-style", action="store_true",
                        help="Load the custom style sheet defined in EXTERNAL_STYLE_PATH.")
    parser.add_argument("--style", type=str, default=None,
                        help="Built-in Matplotlib stylesheet name (e.g., 'seaborn-v0_8-whitegrid'). Ignored if --use-external-style is set.")
    parser.add_argument("--figsize", nargs=2, type=float, default=[8.0, 6.0],
                        metavar=("WIDTH", "HEIGHT"), help="Base figure size in inches per subfigure (e.g., --figsize 10 8).")
    parser.add_argument("--dpi", type=int, default=150, help="DPI resolution for saved figures.")
    parser.add_argument("--list", action="store_true", default=False, help="list available runs")

    args = parser.parse_args()

    available_runs = discover_runs(args.data_dir)
    if not available_runs:
        print(f"Error: No successful run data found in {args.data_dir} (Checked GridFire/Ok/).")
        sys.exit(1)

    if args.list:
        for run in available_runs:
            print(f"==> {run}")
        exit()

    if "all" in [r.lower() for r in args.runs]:
        target_runs = available_runs
    else:
        target_runs = [r for r in args.runs if r in available_runs]
        missing = [r for r in args.runs if r.lower() != "all" and r not in target_runs]
        if missing:
            print(f"Warning: The following runs were skipped because they failed or weren't found: {', '.join(missing)}")

    if args.plot and target_runs:
        fig_opts = {
            "figsize": tuple(args.figsize),
            "dpi": args.dpi,
            "style": args.style,
            "use_ext_style": args.use_external_style,
            "merge_runs": args.merge_runs,
            "x_scale": args.x_scale,
            "y_scale": args.y_scale
        }
        plot_data(target_runs, args.plot, args.data_dir, args.compare_pynucastro, args.format, fig_opts)
    elif args.plot:
        print("Error: No valid runs matched your selection.")

if __name__ == "__main__":
    main()