import os.path

from gridfire.engine import GraphEngine, NetworkJacobian
from fourdst.composition import Composition
from gridfire._gridfire.engine.scratchpads import StateBlob
from gridfire._gridfire.engine.scratchpads import ScratchPadType
import gridfire
import numpy as np

from typing import Tuple
from tqdm import tqdm

def get_stiffness_ratio_and_eigs(J: NetworkJacobian) -> Tuple[float, np.ndarray]:
    jac = J.to_numpy()
    eigenvalues = np.linalg.eigvals(jac)
    abs_real_eigvals = np.abs(np.real(eigenvalues))
    non_zero_eigenvalues = abs_real_eigvals[abs_real_eigvals > 0]
    
    if len(non_zero_eigenvalues) == 0:
        return 1.0, eigenvalues

    min_lambda = np.min(non_zero_eigenvalues)
    max_lambda = np.max(non_zero_eigenvalues)

    if (min_lambda) == 0:
        return np.inf, eigenvalues

    stiffness_ratio = max_lambda / min_lambda
    return stiffness_ratio, eigenvalues

def setup() -> Tuple[Composition, GraphEngine, StateBlob]:
    baseComposition : Composition = Composition({"H-1": 0.702, "He-4": 0.06, "O-16": 1e-5})
    engine: GraphEngine = GraphEngine(baseComposition, 4)
    blob = StateBlob()
    blob.enroll(ScratchPadType.GRAPH_ENGINE_SCRATCHPAD)

    graph_ctx = blob.get(ScratchPadType.GRAPH_ENGINE_SCRATCHPAD)
    graph_ctx.initialize(engine)
    return baseComposition, engine, blob

def format_grid_point_stiff_data(J: NetworkJacobian, engine: GraphEngine, blob: StateBlob, T9: float, density: float) -> str:
    species = engine.getNetworkSpecies(blob)
    result = []
    for sp_row in species:
        for sp_col in species:
            Jij = J[sp_row, sp_col]
            if Jij != 0:
                result.append(f"J[{sp_row},{sp_col}](T9={T9},rho={density}) = {Jij}")
    return "\n".join(result)

def get_stiff(engine: GraphEngine, blob: StateBlob, comp: Composition, T9: float, density: float) -> Tuple[float, np.ndarray]:
    J = engine.generateJacobianMatrix(blob, comp, T9, density)
    S, eigs = get_stiffness_ratio_and_eigs(J)
    return S, eigs

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    comp, engine, blob = setup()

    if os.path.exists("stiff_species.dat"):
        os.remove("stiff_species.dat")

    T = np.linspace(4e6, 1e9, 20)
    R = np.logspace(2, 6, 20)

    TT, RR = np.meshgrid(T, R)

    print(f"Stiffness for a T9=0.015: {get_stiff(engine, blob, comp, 0.015, 160)[0]}")
    print(f"Stiffness for a T9=0.040: {get_stiff(engine, blob, comp, 0.040, 160)[0]}")

    # Generate the global stiffness map
    SS = np.zeros_like(TT)
    for tid, t in tqdm(enumerate(T), total=len(T), desc="Temperature Grid"):
        for rid, r in tqdm(enumerate(R), total=len(R), desc=f"Density Grid at T={t:5.3e}", leave=False):
            stiffness, _ = get_stiff(engine, blob, comp, t/1e9, r)
            SS[tid, rid] = np.log10(stiffness)

    np.savetxt("SS.np.dat", SS)
    np.savetxt("T.np.dat", T)
    np.savetxt("R.np.dat", R)
    with open("metadata.txt", "w") as f:
        f.write(f"SS - Stiffness ratio for jacobian matrix generated from a fully constructed GraphEngine with gridfire version: {gridfire.__version__}\n")
        f.write("T - Temperature [K], first axis (rows)\n")
        f.write("R - Log10 Density [log (g cm^-3)], second axis (cols)\n")
    
    # Generate data specific to the requested T9 histograms
    target_T9s = [0.015, 0.02, 0.1]
    hist_data_reals = []
    valid_mins, valid_maxs = [], []

    for t9 in target_T9s:
        eigs_for_t9 = []
        for r in R:
            _, eigs = get_stiff(engine, blob, comp, t9, r)
            eigs_for_t9.append(eigs)
        
        flat_eigs = np.concatenate(eigs_for_t9)
        real_eigs = np.real(flat_eigs)
        non_zero_reals = real_eigs[real_eigs != 0]
        abs_reals = np.abs(non_zero_reals)
        hist_data_reals.append(abs_reals)
        
        if len(abs_reals) > 0:
            valid_mins.append(np.min(abs_reals))
            valid_maxs.append(np.max(abs_reals))

    # ------------------ Plotting ------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Plot 1: Stiffness Ratio Contour
    img = ax1.contourf(T, R, SS)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    
    cbar = plt.colorbar(img, ax=ax1)
    cbar.ax.invert_yaxis()
    ax1.set_xlabel("Temperature [K]")
    ax1.set_ylabel("Density [g cm$^{-3}$]")
    ax1.set_title("Log10 Stiffness Ratio")

    # Plot 2: Overlaid Log-Spaced Histograms for specific T9 values
    if valid_mins and valid_maxs:
        global_min = np.min(valid_mins)
        global_max = np.max(valid_maxs)
        # Create global bins so the histograms align perfectly
        bins = np.logspace(np.log10(global_min), np.log10(global_max), 50)
        
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c'] # Blue, Orange, Green
        
        for i, t9 in enumerate(target_T9s):
            ax2.hist(hist_data_reals[i], bins=bins, alpha=0.6, label=f"T9 = {t9}", 
                     color=colors[i], edgecolor='black', linewidth=0.5)
            
        ax2.set_xscale('log')
        ax2.legend()
        
    ax2.grid(True, linestyle='--', alpha=0.6, axis='y')
    ax2.set_xlabel(r"Absolute Real Part $|Re(\lambda)|$")
    ax2.set_ylabel("Frequency")
    ax2.set_title("Eigenvalue Distribution at Selected Temperatures")

    plt.tight_layout()
    plt.show()