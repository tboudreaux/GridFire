from enum import Enum

from typing import Dict, List, Any, SupportsFloat
import json
from datetime import datetime
import os
import sys

from gridfire.solver import CVODETimestepContext
import gridfire

class LogEntries(Enum):
    Step = "Step"
    t = "t"
    dt = "dt"
    eps = "eps"
    Composition = "Composition"
    ReactionContributions = "ReactionContributions"


class StepLogger:
    def __init__(self):
        self.num_steps : int = 0
        self.steps : List[Dict[LogEntries, Any]] = []

    def log_step(self, ctx : CVODETimestepContext):
        comp_data: Dict[str, SupportsFloat] = {}
        for species in ctx.engine.getNetworkSpecies():
            sid = ctx.engine.getSpeciesIndex(species)
            comp_data[species.name()] = ctx.state[sid]
        entry : Dict[LogEntries, Any] = {
            LogEntries.Step: ctx.num_steps,
            LogEntries.t: ctx.t,
            LogEntries.dt: ctx.dt,
            LogEntries.Composition: comp_data,
        }
        self.steps.append(entry)
        self.num_steps += 1

    def to_json(self, filename: str, **kwargs):
        serializable_steps : List[Dict[str, Any]] = [
            {
                LogEntries.Step.value: step[LogEntries.Step],
                LogEntries.t.value: step[LogEntries.t],
                LogEntries.dt.value: step[LogEntries.dt],
                LogEntries.Composition.value: step[LogEntries.Composition],
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