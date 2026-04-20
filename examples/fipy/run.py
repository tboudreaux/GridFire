import numpy as np
import fipy as fp
import gridfire as gf
import matplotlib.pyplot as plt
from gridfire.type import NetIn
import fourdst

def main():
    Y_solar = [7.0262E-01, 9.7479E-06, 6.8955E-02, 2.5000E-04, 7.8554E-05, 6.0144E-04, 8.1031E-05, 2.1513E-05]
    S = ["H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"]
    base_comp = fourdst.composition.Composition(S, Y_solar)

    stellar_policy = gf.policy.MainSequencePolicy(base_comp)
    construct = stellar_policy.construct()

    point_solver = gf.solver.PointSolver(construct.engine)
    grid_solver = gf.solver.GridSolver(construct.engine, point_solver)
    solver_ctx = gf.solver.GridSolverContext(construct.scratch_blob)

    solver_ctx.zone_completion_logging = False

    species_list = construct.engine.getNetworkSpecies(construct.scratch_blob)

    plot_species = ["H-1", "He-4", "C-12", "O-16", "N-14", "Mg-24"]
    nx = 100
    dx = 1.0
    mesh = fp.Grid1D(nx=nx, dx=dx)
    center_idx = nx // 2

    T9 = fp.CellVariable(name="Temperature (T9)", mesh=mesh, value=0.015)

    x = mesh.cellCenters[0]
    T9.setValue(0.015 + 0.005 * fp.numerix.exp(-((x - (nx*dx/2))**2) / 50.0))

    species_vars = {}
    for i, sp in enumerate(species_list):
        sp_name = sp.name()
        val = base_comp.getMolarAbundance(sp) if base_comp.contains(sp) else 0.0
        species_vars[sp_name] = fp.CellVariable(name=sp_name, mesh=mesh, value=val)
    D_T = 1e4
    D_Y = 1e2

    eq_T = fp.TransientTerm() == fp.DiffusionTerm(coeff=D_T)
    eqs_Y = {sp_name: fp.TransientTerm() == fp.DiffusionTerm(coeff=D_Y) for sp_name in species_vars}
    dt = 0.1
    t_final = 3e14
    t = 0.0
    rho_const = 160.0
    cp_const = 2.0e8

    plt.ion()
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))

    line_T, = ax1.plot(mesh.cellCenters[0], T9.value, label="T9 (Billion K)", color='red')
    ax1.set_ylabel("Temperature (T9)")
    ax1.legend()

    lines_Y = {}
    colors = plt.cm.tab10(np.linspace(0, 1, len(plot_species)))
    for sp_name, color in zip(plot_species, colors):
        if sp_name in species_vars:
            lines_Y[sp_name], = ax2.plot(mesh.cellCenters[0], species_vars[sp_name].value, label=sp_name, color=color)
    ax2.set_ylabel("Molar Abundance")
    ax2.set_xlabel("Zone Index")
    ax2.set_yscale('log')
    ax2.set_ylim(1e-12, 1.5)
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    time_history = []
    history_vars = {sp: [] for sp in plot_species}
    lines_time = {}

    for sp_name, color in zip(plot_species, colors):
        lines_time[sp_name], = ax3.plot([], [], label=f"Total {sp_name}", color=color, linewidth=2)

    ax3.set_xlabel("Time (s)")
    ax3.set_ylabel("Total Moles (mol/cm$^2$)")
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlim(1e-1, t_final)
    ax3.set_ylim(1e-5, 1e5)
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()

    while t < t_final:
        for sp_name, sp_var in species_vars.items():
            eqs_Y[sp_name].solve(var=sp_var, dt=dt)

        net_ins = []
        T9_array = T9.value

        for i in range(nx):
            net_in = NetIn()
            net_in.temperature = max(1e6, T9_array[i] * 1e9) # Convert T9 back to Kelvin
            net_in.density = rho_const
            net_in.tMax = dt

            local_Y = [max(0.0, species_vars[sp.name()].value[i]) for sp in species_list]

            net_in.composition = fourdst.composition.Composition(
                [sp.name() for sp in species_list], local_Y
            )
            net_ins.append(net_in)

        results = grid_solver.evaluate(solver_ctx, net_ins)

        for i in range(nx):
            for sp in species_list:
                sp_name = sp.name()
                if results[i].composition.contains(sp):
                    species_vars[sp_name].value[i] = results[i].composition.getMolarAbundance(sp)
        time_history.append(max(t, 1e-5))

        for sp_name in plot_species:
            if sp_name in species_vars:
                tot_moles = np.sum(species_vars[sp_name].value * rho_const * dx)
                history_vars[sp_name].append(tot_moles)

        line_T.set_ydata(T9_array)
        for sp_name in plot_species:
            if sp_name in species_vars:
                lines_Y[sp_name].set_ydata(species_vars[sp_name].value)
                lines_time[sp_name].set_data(time_history, history_vars[sp_name])

        ax1.set_ylim(min(T9_array)*0.95, max(T9_array)*1.05)

        valid_mins = []
        valid_maxs = []
        for sp_name in plot_species:
            if len(history_vars[sp_name]) > 0:
                arr = np.array(history_vars[sp_name])
                pos_arr = arr[arr > 0]
                if len(pos_arr) > 0:
                    valid_mins.append(np.min(pos_arr))
                valid_maxs.append(np.max(arr))

        if valid_mins and valid_maxs:
            ax3.set_ylim(min(valid_mins) * 0.5, max(valid_maxs) * 2.0)

        if t > 1e-1:
            ax3.set_xlim(1e-1, max(t * 1.5, 1.0))

        fig.canvas.draw()
        fig.canvas.flush_events()

        t += dt
        print(f"Time: {t:.2e}s | Center Temp: {T9.value[center_idx]*1e9:.2e} K | dt: {dt:.2e}s")
        dt = min(dt * 2, 5e12)

    plt.ioff()
    plt.show()

if __name__ == "__main__":
    main()
