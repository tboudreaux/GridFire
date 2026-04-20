import shutil
from abc import ABC, abstractmethod

import fourdst.atomic
import scipy.integrate
import gridfire
from fourdst.composition import Composition
from gridfire.engine import DynamicEngine, GraphEngine, AdaptiveEngineView, MultiscalePartitioningEngineView
from gridfire.engine import EngineTypes
from gridfire.policy import MainSequencePolicy
from gridfire.type import NetIn, NetOut
from gridfire.exceptions import GridFireError
from gridfire.solver import PointSolver, PointSolverContext
from logger import StepLogger
from typing import List
import re
from typing import Dict, Tuple, Any, Union
from datetime import datetime

import pynucastro as pyna
import os
import importlib.util
import sys
import numpy as np
import json
import time

EngineTypeLookup : Dict[EngineTypes, Any] = {
    EngineTypes.ADAPTIVE_ENGINE_VIEW: AdaptiveEngineView,
    EngineTypes.MULTISCALE_PARTITIONING_ENGINE_VIEW: MultiscalePartitioningEngineView,
    EngineTypes.GRAPH_ENGINE: GraphEngine
}

def load_network_module(filepath):
    module_name = os.path.basename(filepath).replace(".py", "")
    if module_name in sys.modules: # clear any existing module with the same name
        del sys.modules[module_name]

    spec = importlib.util.spec_from_file_location(module_name, filepath)
    if spec is None:
        raise FileNotFoundError(f"Could not find module at {filepath}")

    network_module = importlib.util.module_from_spec(spec)

    sys.modules[module_name] = network_module

    spec.loader.exec_module(network_module)

    return network_module

def get_pyna_rate(my_rate_str, library):
    match = re.match(r"([a-zA-Z0-9]+)\(([^,]+),([^)]*)\)(.*)", my_rate_str)

    if not match:
        print(f"Could not parse string format: {my_rate_str}")
        return None

    target = match.group(1)
    projectile = match.group(2)
    ejectiles = match.group(3)
    product = match.group(4)

    def expand_species(s_str):
        if not s_str or s_str.strip() == "":
            return []

        # Split by space (handling "p a" or "2p a")
        parts = s_str.split()
        expanded = []

        for p in parts:
            # Check for multipliers like 2p, 3a
            mult_match = re.match(r"(\d+)([a-zA-Z0-9]+)", p)
            if mult_match:
                count = int(mult_match.group(1))
                spec = mult_match.group(2)
                # Map common aliases if necessary (though pyna handles most)
                if spec == 'a': spec = 'he4'
                expanded.extend([spec] * count)
            else:
                spec = p
                if spec == 'a': spec = 'he4'
                expanded.append(spec)
        return expanded

    reactants_str = [target] + expand_species(projectile)
    products_str = expand_species(ejectiles) + [product]

    # Convert strings to pyna.Nucleus objects
    try:
        r_nuc = [pyna.Nucleus(r) for r in reactants_str]
        p_nuc = [pyna.Nucleus(p) for p in products_str]
    except Exception as e:
        print(f"Error converting nuclei for {my_rate_str}: {e}")
        return None

    rates = library.get_rate_by_nuclei(r_nuc, p_nuc)

    if rates:
        if isinstance(rates, list):
            return rates[0]  # Return the first match
        return rates
    else:
        return None

class TestSuite(ABC):
    def __init__(self, name: str, description: str, temp: float, density: float, tMax: float, composition: Composition, notes: str = ""):
        self.name : str = name
        self.description : str = description
        self.temperature : float = temp
        self.density : float = density
        self.tMax : float = tMax
        self.composition : Composition = composition
        self.notes : str = notes

    def evolve_pynucastro(self, engine: DynamicEngine, ctx: PointSolverContext, output: str = "pynucastro"):
        print("Evolution complete. Now building equivalent pynucastro network...")
        # Build equivalent pynucastro network for comparison
        reaclib_library : pyna.ReacLibLibrary = pyna.ReacLibLibrary()
        rate_names = [r.id().replace("e+","").replace("e-","").replace(", ", ",") for r in engine.getNetworkReactions(ctx.engine_ctx)]

        with open(f"{self.name}_rate_names_pynuc.txt", "w") as f:
            for r_name in rate_names:
                f.write(f"{r_name}\n")

        goodRates : List[pyna.rates.reaclib_rate.ReacLibRate] = []
        missingRates = []

        for r_str in rate_names:
            # Try the exact name match first (fastest)
            try:
                pyna_rate = reaclib_library.get_rate_by_name(r_str)
                if isinstance(pyna_rate, list):
                    goodRates.append(pyna_rate[0])
                else:
                    goodRates.append(pyna_rate)
            except:
                # Fallback to the smart parser
                pyna_rate = get_pyna_rate(r_str, reaclib_library)
                if pyna_rate:
                    goodRates.append(pyna_rate)
                else:
                    missingRates.append(r_str)


        pynet : pyna.PythonNetwork = pyna.PythonNetwork(rates=goodRates)
        pynet.write_network(f"{self.name}_pynucastro_network.py")

        net = load_network_module(f"{self.name}_pynucastro_network.py")
        Y0 = np.zeros(net.nnuc)
        Y0[net.jp] = self.composition.getMolarAbundance("H-1")
        Y0[net.jhe3] = self.composition.getMolarAbundance("He-3")
        Y0[net.jhe4] = self.composition.getMolarAbundance("He-4")
        Y0[net.jc12] = self.composition.getMolarAbundance("C-12")
        Y0[net.jn14] = self.composition.getMolarAbundance("N-14")
        Y0[net.jo16] = self.composition.getMolarAbundance("O-16")
        Y0[net.jne20] = self.composition.getMolarAbundance("Ne-20")
        Y0[net.jmg24] = self.composition.getMolarAbundance("Mg-24")

        print("Starting pynucastro integration...")
        startTime = time.time()
        sol = scipy.integrate.solve_ivp(
            net.rhs,
            [0, self.tMax],
            Y0,
            args=(self.density, self.temperature),
            method="BDF",
            jac=net.jacobian,
            rtol=1e-5,
            atol=1e-8
        )
        endTime = time.time()
        initial_duration = endTime - startTime
        print("Pynucastro integration complete. Writing results to JSON...")
        print("Running pynucastro a second time to account for any JIT compilation overhead...")
        startTime = time.time()
        sol = scipy.integrate.solve_ivp(
            net.rhs,
            [0, self.tMax],
            Y0,
            args=(self.density, self.temperature),
            method="BDF",
            jac=net.jacobian,
            rtol=1e-5,
            atol=1e-8
        )
        endTime = time.time()
        final_duration = endTime - startTime
        print(f"Pynucastro second integration complete. Initial run time: {initial_duration: .4f} s, Second run time: {final_duration: .4f} s")

        data: List[Dict[str, Union[float, Dict[str, float]]]] = []

        for time_step, t in enumerate(sol.t):
            data.append({"t": t, "Composition": {}})
            for j in range(net.nnuc):
                A = net.A[j]
                Z = net.Z[j]
                species: str
                try:
                    species = fourdst.atomic.az_to_species(A, Z).name()
                except:
                    species = f"SP-A_{A}_Z_{Z}"

                data[-1]["Composition"][species] = sol.y[j, time_step]

        pynucastro_json : Dict[str, Any] = {
            "Metadata": {
                "Name": f"{self.name}_pynucastro",
                "Description": f"pynucastro simulation equivalent to GridFire validation suite: {self.description}",
                "Status": "Success",
                "Notes": self.notes,
                "Temperature": self.temperature,
                "Density": self.density,
                "tMax": self.tMax,
                "RunTime0": initial_duration,
                "RunTime1": final_duration,
                "DateCreated": datetime.now().isoformat(),
                "NumSpecies": net.nnuc
            },
            "Steps": data
        }

        filename: str = f"{self.name}_pynucastro.json"
        filepath: str = os.path.join(output, filename)
        with open(filepath, "w") as f:
            json.dump(pynucastro_json, f, indent=4)

    def evolve(self, engine: DynamicEngine, solver_ctx: PointSolverContext, netIn: NetIn, pynucastro_compare: bool = True, engine_type: EngineTypes | None = None, output: str = "output"):
        solver : PointSolver = PointSolver(engine)

        stepLogger : StepLogger = StepLogger()
        solver_ctx.callback = lambda ctx: stepLogger.log_step(ctx)

        startTime = time.time()
        subdir: str = os.path.join(output, "GridFire")
        os.makedirs(subdir, exist_ok=True)
        try:
            startTime = time.time()
            netOut : NetOut = solver.evaluate(solver_ctx, netIn)
            endTime = time.time()
            filename: str = f"{self.name}_OKAY.json"
            subdir2: str = os.path.join(subdir, "Ok")
            os.makedirs(subdir2, exist_ok=True)
            filepath: str = os.path.join(subdir2, filename)
            stepLogger.to_json(
                filepath,
                Name = f"{self.name}_Success",
                Description=self.description,
                Status="Success",
                Notes=self.notes,
                Temperature=netIn.temperature,
                Density=netIn.density,
                tMax=netIn.tMax,
                FinalEps = netOut.energy,
                FinaldEpsdT = netOut.dEps_dT,
                FinaldEpsdRho = netOut.dEps_dRho,
                ElapsedTime = endTime - startTime,
                NumSpecies = engine.getNetworkSpecies(solver_ctx.engine_ctx).__len__(),
                NumReactions = engine.getNetworkReactions(solver_ctx.engine_ctx).__len__()
            )
        except GridFireError as e:
            endTime = time.time()
            filename : str = f"{self.name}_FAIL.json"
            subdir2 : str = os.path.join(subdir, "Err")
            os.makedirs(subdir2, exist_ok=True)
            filepath = os.path.join(subdir2, filename)
            stepLogger.to_json(
                filepath,
                Name = f"{self.name}_Failure",
                Description=self.description,
                Status=f"Error",
                ErrorMessage=str(e),
                Notes=self.notes,
                Temperature=netIn.temperature,
                Density=netIn.density,
                tMax=netIn.tMax,
                ElapsedTime = endTime - startTime
            )

        if pynucastro_compare:
            pynuc_dir = os.path.join(output, "pynucastro")
            os.makedirs(pynuc_dir, exist_ok=True)
            if engine_type is not None:
                if engine_type == EngineTypes.MULTISCALE_PARTITIONING_ENGINE_VIEW:
                    print("Pynucastro comparison using MultiscalePartitioningEngineView...")
                    graphEngine : GraphEngine = GraphEngine(self.composition, depth=3)
                    multiScaleEngine : MultiscalePartitioningEngineView = MultiscalePartitioningEngineView(graphEngine)
                    self.evolve_pynucastro(multiScaleEngine, solver_ctx, pynuc_dir)
                elif engine_type == EngineTypes.GRAPH_ENGINE:
                    print("Pynucastro comparison using GraphEngine...")
                    graphEngine : GraphEngine = GraphEngine(self.composition, depth=3)
                    self.evolve_pynucastro(graphEngine, solver_ctx, pynuc_dir)
                else:
                    print(f"Pynucastro comparison not implemented for engine type: {engine_type}")


    @abstractmethod
    def __call__(self, pynucastro_compare: bool = False, pync_engine: str = "GraphEngine", output: str = "output"):
        pass
