from fourdst.composition import Composition
from gridfire.type import NetIn
from gridfire.policy import MainSequencePolicy
from gridfire.solver import CVODESolverStrategy
from enum import Enum
from typing import Dict, Union, SupportsFloat
import json
import dicttoxml

def init_composition() -> Composition:
    Y = [7.0262E-01, 9.7479E-06, 6.8955E-02, 2.5000E-04, 7.8554E-05, 6.0144E-04, 8.1031E-05, 2.1513E-05]
    S = ["H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"]
    return Composition(S, Y)

def init_netIn(temp: float, rho: float, time: float, comp: Composition) -> NetIn:
    netIn = NetIn()
    netIn.temperature = temp
    netIn.density = rho
    netIn.tMax = time
    netIn.dt0 = 1e-12
    netIn.composition = comp
    return netIn

class StepData(Enum):
    TIME = 0
    DT = 1
    COMP = 2
    CONTRIB = 3


class StepLogger:
    def __init__(self):
        self.num_steps: int = 0
        self.step_data: Dict[int, Dict[StepData, Union[SupportsFloat, Dict[str, SupportsFloat]]]] = {}

    def log_step(self, context):
        engine = context.engine
        self.step_data[self.num_steps] = {}
        self.step_data[self.num_steps][StepData.TIME] = context.t
        self.step_data[self.num_steps][StepData.DT] = context.dt
        comp_data: Dict[str, SupportsFloat] = {}
        for species in engine.getNetworkSpecies():
            sid = engine.getSpeciesIndex(species)
            comp_data[species.name()] = context.state[sid]
        self.step_data[self.num_steps][StepData.COMP] = comp_data
        self.num_steps += 1

    def to_json (self, filename: str):
        serializable_data = {
            stepNum: {
                StepData.TIME.name: step[StepData.TIME],
                StepData.DT.name: step[StepData.DT],
                StepData.COMP.name: step[StepData.COMP],
            }
            for stepNum, step in self.step_data.items()
        }
        with open(filename, 'w') as f:
            json.dump(serializable_data, f, indent=4)

    def to_xml(self, filename: str):
        serializable_data = {
            stepNum: {
                StepData.TIME.name: step[StepData.TIME],
                StepData.DT.name: step[StepData.DT],
                StepData.COMP.name: step[StepData.COMP],
            }
            for stepNum, step in self.step_data.items()
        }
        xml_data = dicttoxml.dicttoxml(serializable_data, custom_root='StepLog', attr_type=False)
        with open(filename, 'wb') as f:
            f.write(xml_data)

def main(temp: float, rho: float, time: float):
    comp = init_composition()
    netIn = init_netIn(temp, rho, time, comp)

    policy = MainSequencePolicy(comp)
    engine = policy.construct()

    solver = CVODESolverStrategy(engine)

    step_logger = StepLogger()
    solver.set_callback(lambda context: step_logger.log_step(context))

    solver.evaluate(netIn, False)
    step_logger.to_xml("log_data.xml")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Simple python example of GridFire usage")
    parser.add_argument("-t", "--temp", type=float, help="Temperature in K", default=1.5e7)
    parser.add_argument("-r", "--rho", type=float, help="Density in g/cm^3", default=1.5e2)
    parser.add_argument("--tMax", type=float, help="Time in s", default=3.1536 * 1e17)
    args = parser.parse_args()
    main(args.temp, args.rho, args.tMax)
