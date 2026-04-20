import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import pynucastro as pyna
import os
import sys
import importlib.util
import time
import matplotlib.lines as mlines
import re
import json
import argparse

from fourdst.composition import Composition
from gridfire.type import NetIn
from gridfire.engine import GraphEngine
from gridfire.solver import PointSolver, PointSolverContext

from tqdm import tqdm

from fourdst.composition.utils import buildCompositionFromMassFractions



def T9(age):
    return 10.0 / np.sqrt(age)

def get_density(age):
    return 4e-5 * (T9(age) ** 3)

def get_pyna_rate(my_rate_str, library):
    match = re.match(r"([a-zA-Z0-9]+)\(([^,]*),([^)]*)\)(.*)", my_rate_str)

    if not match:
        print(f"Could not parse string format: {my_rate_str}")
        return []

    target = match.group(1)
    projectile = match.group(2)
    ejectiles = match.group(3)
    product = match.group(4)

    def expand_species(s_str):
        if not s_str or s_str.strip() == "":
            return []

        parts = s_str.split()
        expanded = []

        for p in parts:
            if p == 'g':
                continue

            mult_match = re.match(r"^(\d+)([a-zA-Z0-9]+)$", p)
            if mult_match:
                count = int(mult_match.group(1))
                spec = mult_match.group(2)
            else:
                count = 1
                spec = p

            if spec == 'g':
                continue
            if spec == 'a': spec = 'he4'

            expanded.extend([spec] * count)
        return expanded

    reactants_str = [target] + expand_species(projectile)
    products_str = expand_species(ejectiles) + [product]

    try:
        r_nuc = [pyna.Nucleus(r) for r in reactants_str]
        p_nuc = [pyna.Nucleus(p) for p in products_str]
    except Exception as e:
        print(f"Error converting nuclei for {my_rate_str}: {e}")
        return []

    rates = library.get_rate_by_nuclei(r_nuc, p_nuc)

    if rates:
        if not isinstance(rates, list):
            return [rates]
        return rates

    r_nuc_names = sorted([str(n) for n in r_nuc])
    p_nuc_names = sorted([str(n) for n in p_nuc])

    ignore_list = ['e-', 'e+', 'g', 'nu', 'anu']

    matched_rates = []
    for rate in library.get_rates():
        lib_r_names = sorted([str(n) for n in rate.reactants if str(n) not in ignore_list])
        lib_p_names = sorted([str(n) for n in rate.products if str(n) not in ignore_list])

        if r_nuc_names == lib_r_names and p_nuc_names == lib_p_names:
            matched_rates.append(rate)

    return matched_rates


def load_network_module(filepath):
    module_name = os.path.basename(filepath).replace(".py", "")
    if module_name in sys.modules:
        del sys.modules[module_name]

    spec = importlib.util.spec_from_file_location(module_name, filepath)
    if spec is None:
        raise FileNotFoundError(f"Error: could not find module at {filepath}")

    network_module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = network_module
    spec.loader.exec_module(network_module)
    return network_module

def main(args):
    tMax = 3600.0
    h = 0.01
    current_time = 180.0

    XpXn = 7.17
    Xn = 1.0 / (1.0 + XpXn)
    Xp = 1.0 - Xn
    comp: Composition = buildCompositionFromMassFractions(["H-1", "n-1"], [Xp, Xn])

    netIn = NetIn()
    netIn.composition = comp
    netIn.dt0 = 1e-12

    if args.depth is not None:
        print(f"Initializing GridFire GraphEngine with restricted depth = {args.depth}")
        engine = GraphEngine(comp, args.depth)
    else:
        print("Initializing full-depth GridFire GraphEngine (Note: pynucastro may take a long time to run JIT, set NUMBA_DISABLE_JIT=1 as an eviromental variable to disable JIT, this makes per timestep time increase but may still be faster for large networks due to the lack of upfront compilation time)")
        engine = GraphEngine(comp)

    blob = engine.constructStateBlob()
    solver_ctx = PointSolverContext(blob)
    solver_ctx.stdout_logging = False
    solver = PointSolver(engine)

    gf_initial_Y = {}
    for sp in engine.getNetworkSpecies(solver_ctx.engine_ctx):
        if comp.contains(sp):
            gf_initial_Y[sp.name()] = comp.getMolarAbundance(sp)
        else:
            gf_initial_Y[sp.name()] = 0.0

    gf_time = []
    gf_results = {}
    step_conditions = []

    gf_start_time = time.time()

    gf_current_time = current_time
    total_steps = int(np.ceil(np.log(tMax / current_time) / np.log(1 + h)))
    with tqdm(total=total_steps, desc="GridFire BBN", unit="step") as pbar:
        while gf_current_time < tMax:
            current_dt = h * gf_current_time
            next_time = gf_current_time + current_dt

            burn_temp = (T9(gf_current_time) + T9(next_time)) / 2.0 * 1e9
            burn_density = (get_density(gf_current_time) + get_density(next_time)) / 2.0

            netIn.temperature = burn_temp
            netIn.density = burn_density
            netIn.tMax = current_dt

            netOut = solver.evaluate(solver_ctx, netIn)
            netIn.composition = netOut.composition

            pbar.update(1)
            pbar.set_postfix(t=f"{gf_current_time:.2e}", T=f"{burn_temp:.2e}", rho=f"{burn_density:.2e}")


            step_conditions.append({
                "dt": current_dt,
                "T": burn_temp,
                "rho": burn_density,
                "t": gf_current_time
            })

            gf_time.append(gf_current_time)

            for sp in engine.getNetworkSpecies(solver_ctx.engine_ctx):
                name = sp.name()
                if name not in gf_results:
                    gf_results[name] = []
                gf_results[name].append(netOut.composition.getMolarAbundance(sp))

            gf_current_time += current_dt

    gf_end_time = time.time()
    print(f"GridFire integration finished in {gf_end_time - gf_start_time:.4f} seconds.")


    print("Building Pynucastro BBN Network...")
    reaclib_library = pyna.ReacLibLibrary()

    rate_names = [r.id().replace("e+","").replace("e-","").replace(", ", ",") for r in engine.getNetworkReactions(solver_ctx.engine_ctx)]

    goodRates = []
    missingRates = []
    skipped_photo_rates = 0
    pyna_rate_mapping = {}
    import io
    import contextlib

    for r_str in rate_names:
        pyna_rates_for_reaction = []

        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            try:
                res = reaclib_library.get_rate_by_name(r_str)
                if res is not None:
                    if isinstance(res, list):
                        pyna_rates_for_reaction.extend(res)
                    else:
                        pyna_rates_for_reaction.append(res)
            except:
                pass

        if not pyna_rates_for_reaction:
            res_nuc = get_pyna_rate(r_str, reaclib_library)
            if res_nuc:
                if isinstance(res_nuc, list):
                    pyna_rates_for_reaction.extend(res_nuc)
                else:
                    pyna_rates_for_reaction.append(res_nuc)

        if pyna_rates_for_reaction:
            pyna_rate_mapping[r_str] = pyna_rates_for_reaction
            for rate in pyna_rates_for_reaction:
                if args.filter_photo:
                    is_photo_rate = any(str(r).lower() in ['g', 'gamma'] for r in rate.reactants)
                    if is_photo_rate:
                        skipped_photo_rates += 1
                        continue

                goodRates.append(rate)
        else:
            missingRates.append(r_str)

    if missingRates:
        print(f"Warning: Could not map {len(missingRates)} rates to Pynucastro (likely absent from default ReacLib).")
        print(f"Missing sample: {missingRates[:10]}...")

    if args.filter_photo:
        print(f"Info: Skipped {skipped_photo_rates} photodisintegration rates due to --filter-photo flag.")

    print("--- Evaluating reaction rates over all temperatures ---")
    gf_rates_history = {}
    py_rates_history = {}
    gf_rate_labels = {}
    py_rate_labels = {}

    for reaction in engine.getNetworkReactions(solver_ctx.engine_ctx):
        r_str = reaction.id().replace("e+","").replace("e-","").replace(", ", ",")
        gf_rates_history[r_str] = []
        py_rates_history[r_str] = []

        try:
            gf_rate_labels[r_str] = reaction.sources()
        except AttributeError:
            try:
                gf_rate_labels[r_str] = reaction.sourceLabel()
            except AttributeError:
                gf_rate_labels[r_str] = "Unknown"

        if r_str in pyna_rate_mapping:
            py_rate_labels[r_str] = [getattr(pr, 'label', 'Unknown') for pr in pyna_rate_mapping[r_str]]
        else:
            py_rate_labels[r_str] = []

    for step in tqdm(step_conditions, desc="Calculating Rates", unit="step"):
        T9_val = step["T"] / 1e9
        T_K = step["T"]

        for reaction in engine.getNetworkReactions(solver_ctx.engine_ctx):
            r_str = reaction.id().replace("e+","").replace("e-","").replace(", ", ",")

            gf_rate_val = 0.0
            try:
                gf_rate_val = reaction.calculate_rate(T9_val, 0, [])
            except:
                try:
                    gf_rate_val = reaction.calculate_rate(T9_val, 0, 0, 0, [], dict())
                except Exception as e:
                    pass
            gf_rates_history[r_str].append(gf_rate_val)

            py_rate_val = 0.0
            if r_str in pyna_rate_mapping:
                for pr in pyna_rate_mapping[r_str]:
                    py_rate_val += pr.eval(T_K)
            py_rates_history[r_str].append(py_rate_val)

    print("--- Rate Comparison Summary ---")
    threshold = 1e-4
    mismatches = {}
    for r_str in gf_rates_history:
        gf_arr = np.array(gf_rates_history[r_str])
        py_arr = np.array(py_rates_history[r_str])

        with np.errstate(divide='ignore', invalid='ignore'):
            denom = np.where(py_arr != 0, py_arr, gf_arr)
            denom = np.where(denom == 0, 1e-30, denom)
            rel_diffs = np.abs(gf_arr - py_arr) / denom

        max_diff = np.max(rel_diffs)
        if max_diff > threshold:
            max_idx = np.argmax(rel_diffs)
            mismatches[r_str] = {
                "max_diff": max_diff,
                "temp": step_conditions[max_idx]["T"],
                "gf_val": gf_arr[max_idx],
                "py_val": py_arr[max_idx]
            }

    if mismatches:
        print(f"Found {len(mismatches)} rates with differences > {threshold:.2%}")
        for r_str, info in mismatches.items():
            gf_lbl = gf_rate_labels.get(r_str, 'Unknown')
            py_lbl = py_rate_labels.get(r_str, [])
            print(f"{r_str:20}: Max Diff = {info['max_diff']:.2%}, at T = {info['temp']:.2e} K")
            print(f"                      GF = {info['gf_val']:.4e} (Source: {gf_lbl})")
            print(f"                      Py = {info['py_val']:.4e} (Sources: {py_lbl})")
    else:
        print(f"All rates match within the {threshold:.2%} threshold across all temperatures.")
    print("-------------------------------")

    pynet = pyna.PythonNetwork(rates=goodRates)
    network_file = "pynuc_bbn_network.py"
    pynet.write_network(network_file)

    net = load_network_module(network_file)

    mapping = {
        "H-1": ("p", "tab:blue"),
        "n-1": ("n", "tab:orange"),
        "He-4": ("he4", "tab:green"),
        "H-2": ("d", "tab:red"),
        "H-3": ("t", "tab:purple"),
        "He-3": ("he3", "tab:brown"),
        "Li-7": ("li7", "tab:pink"),
        "Be-7": ("be7", "tab:gray")
    }

    Y0 = np.zeros(net.nnuc)

    for i, nuc in enumerate(pynet.get_nuclei()):
        nuc_name = str(nuc)
        gf_name = None
        for gf, (py, _) in mapping.items():
            if py == nuc_name:
                gf_name = gf
                break

        if not gf_name:
            match = re.match(r"([a-zA-Z]+)(\d+)", nuc_name)
            if match:
                gf_name = f"{match.group(1).capitalize()}-{match.group(2)}"

        if gf_name and gf_name in gf_initial_Y:
            Y0[i] = gf_initial_Y[gf_name]

    pyna_time = []

    pyna_nuc_names = [str(n) for n in pynet.get_nuclei()]
    pyna_results = {nuc: [] for nuc in pyna_nuc_names}

    pyna_start_time = time.time()

    for step in tqdm(step_conditions, unit="step", desc="pynucastro Integration"):
        sol = scipy.integrate.solve_ivp(
            net.rhs,
            [0, step["dt"]],
            Y0,
            args=(step["rho"], step["T"]),
            method="Radau",
            jac=net.jacobian,
            rtol=1e-8,
            atol=1e-20
        )

        Y0 = sol.y[:, -1]
        pyna_time.append(step["t"])

        for j in range(net.nnuc):
            nuc_name = str(pynet.get_nuclei()[j])
            if nuc_name in pyna_results:
                pyna_results[nuc_name].append(Y0[j])

    pyna_end_time = time.time()
    print(f"Pynucastro integration finished in {pyna_end_time - pyna_start_time:.4f} seconds.")


    export_data = {
        "metadata": {
            "tMax": tMax,
            "h": h,
            "initial_time": current_time,
            "initial_XpXn_ratio": XpXn,
            "initial_mass_fractions": {
                "Xp": Xp,
                "Xn": Xn
            },
            "execution_times_seconds": {
                "gridfire": gf_end_time - gf_start_time,
                "pynucastro": pyna_end_time - pyna_start_time
            },
            "missing_pynucastro_rates": missingRates,
            "skipped_photodisintegration_rates": skipped_photo_rates if args.filter_photo else 0,
            "rate_labels": {
                "gridfire": gf_rate_labels,
                "pynucastro": py_rate_labels
            }
        },
        "thermodynamic_conditions": step_conditions,
        "data": {
            "gridfire": {
                "time": gf_time,
                "molar_abundances": gf_results,
                "reaction_rates": gf_rates_history
            },
            "pynucastro": {
                "time": pyna_time,
                "molar_abundances": pyna_results,
                "reaction_rates": py_rates_history
            }
        }
    }

    json_out_file = "bbn_simulation_data.json"
    with open(json_out_file, "w") as f:
        json.dump(export_data, f, indent=4)
    plt.style.use("default")
    fig, ax = plt.subplots(figsize=(10, 7))

    for gf_name, (pyna_name, color) in mapping.items():
        if gf_name in gf_results:
            ax.plot(gf_time, gf_results[gf_name], color=color, linestyle="-", linewidth=2.5, label=f"GF {gf_name}")
        if pyna_name in pyna_results:
            ax.plot(pyna_time, pyna_results[pyna_name], color=color, linestyle="--", linewidth=1.5, label=f"Pyna {pyna_name}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-12, 2)
    ax.set_xlabel("Time (s)", fontsize=14)
    ax.set_ylabel("Molar Abundance (Y)", fontsize=14)

    line_gf = mlines.Line2D([], [], color='black', linestyle='-', linewidth=2.5, label='GridFire')
    line_py = mlines.Line2D([], [], color='black', linestyle='--', linewidth=1.5, label='Pynucastro')

    sp_handles = []
    for gf_name, (pyna_name, color) in mapping.items():
        sp_handles.append(mlines.Line2D([], [], color=color, linestyle='-', linewidth=2, label=gf_name))

    ax.legend(handles=[line_gf, line_py] + sp_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=12)

    out_file = "bbn_comparison.pdf"
    plt.savefig(out_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GridFire vs Pynucastro BBN Comparison")
    parser.add_argument("--filter-photo", action="store_true",
                        help="Filter out photodisintegration (reverse) rates to mimic GridFire's forward-only mechanics.")
    parser.add_argument("--depth", type=int, default=None,
                        help="Limit the assembly depth of GridFire's GraphEngine. E.g., setting '--depth 3' shrinks the network size from 5000+ reactions to ~100, which reduces Pynucastro's Numba JIT compile time from hours to seconds.")
    args = parser.parse_args()

    main(args)
