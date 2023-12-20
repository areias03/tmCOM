from typing import List
import seaborn as sns
import matplotlib as plt

import pandas as pd
import numpy as np
from cobra.io import read_sbml_model
from cobra.io.sbml import Model
from mewpy.simulation import Environment, Simulator, get_simulator

from concurrent.futures import ProcessPoolExecutor


ec_bt = read_sbml_model('../models/ec/ec_Bacteroides_thetaiotaomicron_VPI_5482.xml')
ec_bu = read_sbml_model('../models/ec/ec_Bacteroides_uniformis_ATCC_8492.xml')
ec_ec = read_sbml_model('../models/ec/ec_Escherichia_coli_ED1a.xml')
ec_fn = read_sbml_model('../models/ec/ec_Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586.xml')
ec_ri = read_sbml_model('../models/ec/ec_Roseburia_intestinalis_L1_82.xml')
ec_sp = read_sbml_model('../models/ec/ec_Streptococcus_parasanguinis_ATCC_15912.xml')
ec_ss = read_sbml_model('../models/ec/ec_Streptococcus_salivarius_DSM_20560.xml')
ec_cc = read_sbml_model('../models/ec/ec_Coprococcus_comes_ATCC_27758.xml')

m_list = [ec_bt, ec_bu, ec_cc, ec_ec, ec_fn, ec_ri, ec_sp, ec_ss]


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
    return df


if __name__ == "__main__":

    merged: List = []

    with ProcessPoolExecutor() as p:
        solution_list = []
        for result in p.map(prot_pool_analysis, m_list):
            merged.append(result)
        df = pd.concat(merged, axis=1)
        print(df)
        df_norm = (df - df.mean()) / df.std()
        sns.heatmap(df_norm)
