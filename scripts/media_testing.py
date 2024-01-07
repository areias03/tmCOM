import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from typing import List

from cobra.io import read_sbml_model
from cobra.io.sbml import Model
from mewpy.simulation import Simulator, get_simulator

# Models

ec_bt = read_sbml_model('../models/ec/ec_Bacteroides_thetaiotaomicron_VPI_5482.xml')
ec_bu = read_sbml_model('../models/ec/ec_Bacteroides_uniformis_ATCC_8492.xml')
ec_ec = read_sbml_model('../models/ec/ec_Escherichia_coli_ED1a.xml')
ec_fn = read_sbml_model('../models/ec/ec_Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586.xml')
ec_ri = read_sbml_model('../models/ec/ec_Roseburia_intestinalis_L1_82.xml')
ec_sp = read_sbml_model('../models/ec/ec_Streptococcus_parasanguinis_ATCC_15912.xml')
ec_ss = read_sbml_model('../models/ec/ec_Streptococcus_salivarius_DSM_20560.xml')
ec_cc = read_sbml_model('../models/ec/ec_Coprococcus_comes_ATCC_27758.xml')

model_l = [ec_bt, ec_bu, ec_ec, ec_fn, ec_ri, ec_sp, ec_ss, ec_cc]

# Medias

M1 = {'R_EX_fe2[e]': (-0.033, 1000), 'R_EX_pheme[e]': (0, 1000), 'R_prot_pool_exchange': (0, 0.7)}
M2 = {'R_EX_fe2[e]': (0, 1000), 'R_EX_pheme[e]': (0, 1000), 'R_prot_pool_exchange': (0, 0.7)}
M3 = {'R_EX_fe2[e]': (-0.0082, 1000), 'R_EX_pheme[e]': (-0.0007, 1000), 'R_prot_pool_exchange': (0, 0.7)}
M4 = {'R_EX_fe2[e]': (-0.033, 1000), 'R_EX_pheme[e]': (-0.0027, 1000), 'R_prot_pool_exchange': (0, 0.7)}
M5 = {'R_EX_fe2[e]': (0, 1000), 'R_EX_pheme[e]': (-0.0027, 1000), 'R_prot_pool_exchange': (0, 0.7)}
M6 = {'R_EX_fe2[e]': (-0.125, 1000), 'R_EX_pheme[e]': (0, 1000), 'R_prot_pool_exchange': (0, 0.7)}
M7 = {'R_EX_fe2[e]': (0, 1000), 'R_EX_pheme[e]': (-0.125, 1000), 'R_prot_pool_exchange': (0, 0.7)}
M8 = {'R_EX_fe2[e]': (-2.0, 1000), 'R_EX_pheme[e]': (-0.125, 1000), 'R_prot_pool_exchange': (0, 0.7)}
M9 = {'R_EX_fe2[e]': (-0.125, 1000), 'R_EX_pheme[e]': (-0.125, 1000), 'R_prot_pool_exchange': (0, 0.7)}



def media_growth_test(model: Model):
    media_l = [M1, M2, M3, M4, M5, M6, M7, M8, M9]
    sim: Simulator = get_simulator(model)
    sol_list: List = []
    model_id = str(model.id.strip('M_'))

    for m in media_l:
        res = sim.simulate(constraints=m)
        print(res.objective_value)
        sol_list.append(res.objective_value)

    data = {'Constraints': media_l, model_id: sol_list}
    df = pd.DataFrame(data)
    df = df.set_index('Constraints')
    return df


if __name__ == "__main__":
    merged: List = []

    with ProcessPoolExecutor() as p:
        solution_list = []
        for result in p.map(media_growth_test, model_l):
            merged.append(result)
    df = pd.concat(merged, axis=1)
    df.to_csv('../data/model_tuning/media_growth_test.csv')
    print(df)
