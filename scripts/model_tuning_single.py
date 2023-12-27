from typing import List
import seaborn as sns
import matplotlib as plt

import pandas as pd
import numpy as np
from cobra.io import read_sbml_model
from cobra.io.sbml import Model
from mewpy.simulation import Environment, Simulator, get_simulator

import argparse

parser = argparse.ArgumentParser(description='Tuner for ec models.')

parser.add_argument('filepath', type=str, help='Path to a EC model')
parser.add_argument('outpath', type=str, help='Output path. Defaults to ../data/model_tuning/', nargs="?", default = '../data/model_tuning/')

args = parser.parse_args()

model_path = args.filepath
out_path = args.outpath

model = read_sbml_model(model_path)

def prot_pool_analysis(model: Model):
    l_va = np.linspace(0.1, 1.0, num=10)
    constraint_list: List = []
    sol_list: List = []
    ind_list: List = []

    for i in range(len(l_va)):
        constraint_list.append({'prot_pool_exchange': (0, round(l_va[i], 1))})
        ind_list.append(f'0 - {round(l_va[i], 1)}')

    model_id = str(model.id.strip('M_'))
    sim: Simulator = get_simulator(model)
    env = Environment.complete(sim, max_uptake=1000.0, inplace=False)
    sim.set_environmental_conditions(env)
    for i in range(len(l_va)):
        cons = constraint_list[i]
        solution = sim.simulate(constraints=cons)
        sol_list.append(solution.objective_value)
    data = {'Protein pool exchange': ind_list, model_id: sol_list}
    df = pd.DataFrame(data)
    df = df.set_index('Protein pool exchange')
    df.to_csv(str(out_path + model_id)+'.csv')
    return df


if __name__ == "__main__":
    print(prot_pool_analysis(model))
