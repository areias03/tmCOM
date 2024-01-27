import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from typing import List

from cobra.io import read_sbml_model
from cobra.io.sbml import Model
from mewpy.simulation import Simulator, get_simulator

# EC models

ec_bt = read_sbml_model('../models/ec/ec_Bacteroides_thetaiotaomicron_VPI_5482.xml')
ec_bu = read_sbml_model('../models/ec/ec_Bacteroides_uniformis_ATCC_8492.xml')
ec_ec = read_sbml_model('../models/ec/ec_Escherichia_coli_ED1a.xml')
ec_fn = read_sbml_model('../models/ec/ec_Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586.xml')
ec_ri = read_sbml_model('../models/ec/ec_Roseburia_intestinalis_L1_82.xml')
ec_sp = read_sbml_model('../models/ec/ec_Streptococcus_parasanguinis_ATCC_15912.xml')
ec_ss = read_sbml_model('../models/ec/ec_Streptococcus_salivarius_DSM_20560.xml')


# Medias

M1 = {'EX_fe2(e)_REV': (0, 0.000000033)}
M2 = {'EX_fe2(e)_REV': (0, 0), 'EX_pheme(e)_REV': (0, 0),
      'DM_pheme(c)': (0, 0)}
M3 = {'EX_fe2(e)_REV': (0, 0.0000000082), 'EX_pheme(e)_REV': (
    0, 0.0000000007), 'DM_pheme(c)': (0, 0.0000000007)}
M4 = {'EX_fe2(e)_REV': (0, 0.000000033), 'EX_pheme(e)_REV': (
    0, 0.0000000027), 'DM_pheme(c)': (0, 0.0000000027)}
M5 = {'EX_pheme(e)_REV': (0, 0.0000000027), 'DM_pheme(c)': (
    0, 0.0000000027)}
M6 = {'EX_fe2(e)_REV': (0, 0.000000125)}
M7 = {'EX_pheme(e)_REV': (0, 0.000000125), 'DM_pheme(c)': (
    0, 0.000000125)}
M8 = {'EX_fe2(e)_REV': (0, 0.000002), 'EX_pheme(e)_REV': (
    0, 0.000000125), 'DM_pheme(c)': (0, 0.000000125)}
M9 = {'EX_fe2(e)_REV': (0, 0.000000125), 'EX_pheme(e)_REV': (
    0, 0.000000125), 'DM_pheme(c)': (0, 0.000000125)}
M10 = {}

prot_pool = {'M_Bacteroides_thetaiotaomicron_VPI_5482': (0,0.3),
              'M_Bacteroides_uniformis_ATCC_8492': (0,0.3),
              'M_Escherichia_coli_ED1a': (0,0.13),
              'M_Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586': (0,0.32),
              'M_Roseburia_intestinalis_L1_82': (0,0.26),
              'M_Streptococcus_parasanguinis_ATCC_15912': (0,0.2),
              'M_Streptococcus_salivarius_DSM_20560': (0,0.3)}


def find_iron_reactions(model: Model):
    fe_reaction = 'EX_fe2(e)_REV'
    heme_exchange = 'EX_pheme(e)_REV'
    heme_demand = 'DM_pheme(c)'
    has_heme_exchange = False
    has_fe_exchange = False
    has_heme_demand = False
    media_l = [M1, M2, M3, M4, M5, M6, M7, M8, M9, M10]
    for r in model.reactions:
        if 'exchange reaction for heme' in r.name:
            has_heme_exchange = True
        elif 'Exchange of Iron (Fe2+)' in r.name:
            has_fe_exchange = True
        elif 'Demand reaction for heme' in r.name:
            has_heme_demand = True
    for m in media_l:
        if not has_fe_exchange:
            if fe_reaction in m.keys():
                del m[fe_reaction]
        if not has_heme_exchange:
            if heme_exchange in m.keys():
                del m[heme_exchange]
        if not has_heme_demand:
            if heme_demand in m.keys():
                del m[heme_demand]
            


def media_growth_test(model: Model):
    media_name_l = ['33 uM Fe', '0 uM Fe/H',
                    '0,7 uM H; 8,2 uM Fe', '2,7 uM H; 33 uM Fe',
                    '2,7 uM H', '125 uM Fe',
                    '125 uM H', '125 uM H; 2 mM Fe',
                    '125 uM H; 125 uM Fe', 'No constraints']
    sim: Simulator = get_simulator(model)
    sol_list: List = []
    model_id = str(model.id.strip('M_'))
    media_l = [M1, M2, M3, M4, M5, M6, M7, M8, M9, M10]
    for m in media_l:
        m['prot_pool_exchange'] = prot_pool[model.id]
        print(model_id,'\t', m)
        res = sim.simulate(constraints=m)
        sol_list.append(res.objective_value)

    data = {'Constraints': media_name_l, model_id: sol_list}
    df = pd.DataFrame(data)
    df = df.set_index('Constraints')
    print(df)
    return df


def main(model: Model):
    find_iron_reactions(model)
    df = media_growth_test(model)
    return df


if __name__ == "__main__":
    merged: List = []
    model_l = [bt, bu, ec, fn, ri, sp, ss]
    with ProcessPoolExecutor() as p:
        solution_list = []
        for result in p.map(main, model_l):
            merged.append(result)
    df = pd.concat(merged, axis=1)
    df.to_csv('../data/media_testing/ec_media_growth_test.csv')
    print(df)
