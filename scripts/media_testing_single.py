import pandas as pd
from typing import List

from cobra.io import read_sbml_model
from cobra.io.sbml import Model
from mewpy.simulation import Simulator, get_simulator

import argparse

parser = argparse.ArgumentParser(description='Tuner for ec models.')

parser.add_argument('filepath', type=str, help='Path to a EC model')
parser.add_argument('outpath', type=str, help='Output path. Defaults to ../data/media_testing/', nargs="?", default='../data/media_testing/')

args = parser.parse_args()

model_path = args.filepath
out_path = args.outpath

model = read_sbml_model(model_path)

# Medias

M1 = {'EX_fe2(e)_REV': (0, 0.000000033), 'prot_pool_exchange': (0, 1.0)}
M2 = {'EX_fe2(e)_REV': (0, 0), 'EX_pheme(e)_REV': (0, 0), 'DM_pheme(c)': (0, 0), 'prot_pool_exchange': (0, 1.0)}
M3 = {'EX_fe2(e)_REV': (0, 0.0000000082), 'EX_pheme(e)_REV': (0, 0.0000000007), 'DM_pheme(c)': (0, 0.0000000007), 'prot_pool_exchange': (0, 1.0)}
M4 = {'EX_fe2(e)_REV': (0, 0.000000033), 'EX_pheme(e)_REV': (0, 0.0000000027), 'DM_pheme(c)': (0, 0.0000000027), 'prot_pool_exchange': (0, 1.0)}
M5 = {'EX_pheme(e)_REV': (0, 0.0000000027), 'DM_pheme(c)': (0, 0.0000000027), 'prot_pool_exchange': (0, 1.0)}
M6 = {'EX_fe2(e)_REV': (0, 0.000000125), 'prot_pool_exchange': (0, 1.0)}
M7 = {'EX_pheme(e)_REV': (0, 0.000000125), 'DM_pheme(c)': (0, 0.000000125), 'prot_pool_exchange': (0, 1.0)}
M8 = {'EX_fe2(e)_REV': (0, 0.000002), 'EX_pheme(e)_REV': (0, 0.000000125), 'DM_pheme(c)': (0, 0.000000125), 'prot_pool_exchange': (0, 1.0)}
M9 = {'EX_fe2(e)_REV': (0, 0.000000125), 'EX_pheme(e)_REV': (0, 0.000000125), 'DM_pheme(c)': (0, 0.000000125), 'prot_pool_exchange': (0, 1.0)}
M10 = {'prot_pool_exchange': (0, 1.0)}


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
    media_l = [M1, M2, M3, M4, M5, M6, M7, M8, M9, M10]
    sol_list: List = []
    model_id = str(model.id.strip('M_'))
    for m in media_l:
        sim: Simulator = get_simulator(model)
        # sim.set_environmental_conditions(m)
        res = sim.simulate(constraints=m)
        print(res.find(['fe2', 'pheme', 'Iron']))
        print(res.objective_value)
        sol_list.append(res.objective_value)
    data = {'Constraints': media_l, model_id: sol_list}
    df = pd.DataFrame(data)
    df = df.set_index('Constraints')
    df.to_csv(str(out_path + model_id)+'.csv')
    return df


if __name__ == "__main__":
    # for r in model.reactions:
    #     if any(x in r.name for x in ['heme', 'iron']):
    #         print(f'Id: {r.id}\tName: {r.name}')
    find_iron_reactions(model)
    print(media_growth_test(model))
