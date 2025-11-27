from abc import ABC, abstractmethod

import fourdst.atomic
import scipy.integrate
from fourdst.composition import Composition
from gridfire.engine import DynamicEngine, GraphEngine
from gridfire.type import NetIn, NetOut
from gridfire.exceptions import GridFireError
from gridfire.solver import CVODESolverStrategy
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

    def evolve_pynucastro(self, engine: GraphEngine):
        print("Evolution complete. Now building equivalent pynucastro network...")
        # Build equivalent pynucastro network for comparison
        reaclib_library : pyna.ReacLibLibrary = pyna.ReacLibLibrary()
        rate_names = [r.id().replace("e+","").replace("e-","").replace(", ", ",") for r in engine.getNetworkReactions()]

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
        print("Pynucastro integration complete. Writing results to JSON...")

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
                "ElapsedTime": endTime - startTime,
                "DateCreated": datetime.now().isoformat()
            },
            "Steps": data
        }

        with open(f"GridFireValidationSuite_{self.name}_pynucastro.json", "w") as f:
            json.dump(pynucastro_json, f, indent=4)

    def evolve(self, engine: GraphEngine, netIn: NetIn, pynucastro_compare: bool = True):
        solver : CVODESolverStrategy = CVODESolverStrategy(engine)

        stepLogger : StepLogger = StepLogger()
        solver.set_callback(lambda ctx: stepLogger.log_step(ctx))

        startTime = time.time()
        try:
            startTime = time.time()
            netOut : NetOut = solver.evaluate(netIn)
            endTime = time.time()
            stepLogger.to_json(
                f"GridFireValidationSuite_{self.name}_OKAY.json",
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
                ElapsedTime = endTime - startTime
            )
        except GridFireError as e:
            endTime = time.time()
            stepLogger.to_json(
                f"GridFireValidationSuite_{self.name}_FAIL.json",
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
            self.evolve_pynucastro(engine)



    @abstractmethod
    def __call__(self):
        pass
