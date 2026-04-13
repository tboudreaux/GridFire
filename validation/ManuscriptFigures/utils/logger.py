from enum import Enum

from typing import Dict, List, Any, SupportsFloat
import json
from datetime import datetime
import os
import sys

from gridfire.solver import PointSolverTimestepContext
from gridfire._gridfire.engine.scratchpads import StateBlob

from fourdst.composition import Composition
import gridfire

import numpy as np
import pandas as pd

class LogEntries(Enum):
    Step = "Step"
    t = "t"
    dt = "dt"
    eps = "eps"
    Composition = "Composition"
    ReactionContributions = "ReactionContributions"
    MassFractions = "MassFractions"


class StepLogger:
    def __init__(self):
        self.num_steps : int = 0
        self.steps : List[Dict[LogEntries, Any]] = []

    # def log_step(self, ctx: PointSolverTimestepContext):
    #     comp_data: Dict[str, SupportsFloat] = {}
    #     for species in ctx.engine.getNetworkSpecies(ctx.state_ctx):
    #         sid = ctx.engine.getSpeciesIndex(ctx.state_ctx, species)
    #         comp_data[species.name()] = ctx.state[sid]
    #     entry : Dict[LogEntries, Any] = {
    #         LogEntries.Step: ctx.num_steps,
    #         LogEntries.t: ctx.t,
    #         LogEntries.dt: ctx.dt,
    #         LogEntries.eps: ctx.state[-1],
    #         LogEntries.Composition: comp_data,
    #     }
    #     self.steps.append(entry)
    #     self.num_steps += 1

    def log_step(self, ctx: PointSolverTimestepContext):
        full_comp = ctx.composition

        comp_data: Dict[str, float] = {}
        mass_frac: Dict[str, float] = {}
        for species in full_comp.getRegisteredSpecies():
            comp_data[species.name()] = full_comp.getMolarAbundance(species)
            mass_frac[species.name()] = full_comp.getMassFraction(species)

        rhs_calc = ctx.engine.getMostRecentRHSCalculation(ctx.state_ctx)
        instantaneous_eps = rhs_calc.energy if rhs_calc else 0.0

        entry : Dict[LogEntries, Any] = {
            LogEntries.Step: ctx.num_steps,
            LogEntries.t: ctx.t,
            LogEntries.dt: ctx.dt,
            LogEntries.eps: instantaneous_eps,
            LogEntries.Composition: comp_data,
            LogEntries.MassFractions: mass_frac,
        }
        self.steps.append(entry)
        self.num_steps += 1

    def to_json(self, filename: str, **kwargs):
        serializable_steps : List[Dict[str, Any]] = [
            {
                LogEntries.Step.value: step[LogEntries.Step],
                LogEntries.t.value: step[LogEntries.t],
                LogEntries.dt.value: step[LogEntries.dt],
                LogEntries.eps.value: step[LogEntries.eps],
                LogEntries.Composition.value: step[LogEntries.Composition],
                LogEntries.MassFractions.value: step[LogEntries.MassFractions],
            }
            for step in self.steps
        ]
        out_data : Dict[str, Any] = {
            "Metadata": {
                "NumSteps": self.num_steps,
                **kwargs,
                "DateCreated": datetime.now().isoformat(),
                "GridFireVersion": gridfire.__version__,
                "Author": "Emily M. Boudreaux",
                "OS": os.uname().sysname,
                "ClangVersion": os.popen("clang --version").read().strip(),
                "GccVersion": os.popen("gcc --version").read().strip(),
                "PythonVersion": sys.version,
            },
            "Steps": serializable_steps
        }
        with open(filename, 'w') as f:
            json.dump(out_data, f, indent=4)

    @property
    def t(self) -> np.ndarray:
        return np.array([step[LogEntries.t] for step in self.steps])

    @property
    def df(self) -> pd.DataFrame:
        if not self.steps:
            return pd.DataFrame()

        flat_data = []
        for step in self.steps:
            row = {
                "Step": step[LogEntries.Step],
                "t": step[LogEntries.t],
                "dt": step[LogEntries.dt],
                "eps": step[LogEntries.eps],
            }
            X_dict = {f"X_{sp}": x for sp, x in step[LogEntries.MassFractions].items()}
            row.update(step[LogEntries.Composition])
            row.update(X_dict)
            flat_data.append(row)

        df = pd.DataFrame(flat_data)

        df = df.ffill().fillna(0.0)
        return df

    def summary(self) -> Dict[str, Any]:
        if not self.steps:
            return {}
        final_step = self.steps[-1]
        summary_data : Dict[str, Any] = {
            "TotalSteps": self.num_steps,
            "FinalTime": final_step[LogEntries.t],
            "FinalComposition": final_step[LogEntries.Composition],
        }
        return summary_data

    def reset(self):
        self.num_steps = 0
        self.steps : List[Dict[LogEntries, Any]] = []

