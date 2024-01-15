from cobra.io import read_sbml_model
from mewpy.cobra.com import Simulator, SteadyCom, SteadyComVA, mu_score, mp_score, mip_score, mro_score, sc_score, minimal_medium
from mewpy.model.com import CommunityModel
from mewpy.solvers.solution import SimulationResult

import numpy as np
import pandas as pd

from typing import List
from concurrent.futures import ProcessPoolExecutor

# Non-EC models

bt = read_sbml_model('../models/non-ec/agora/Bacteroides_thetaiotaomicron_VPI_5482.xml')
bu = read_sbml_model('../models/non-ec/agora/Bacteroides_uniformis_ATCC_8492.xml')
# ec = read_sbml_model('../models/non-ec/agora/Escherichia_coli_ED1a.xml')
# fn = read_sbml_model('../models/non-ec/agora/Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586.xml')
# ri = read_sbml_model('../models/non-ec/agora/Roseburia_intestinalis_L1_82.xml')
# sp = read_sbml_model('../models/non-ec/agora/Streptococcus_parasanguinis_ATCC_15912.xml')
# ss = read_sbml_model('../models/non-ec/agora/Streptococcus_salivarius_DSM_20560.xml')

# EC models

ec_bt = read_sbml_model('../models/ec/ec_Bacteroides_thetaiotaomicron_VPI_5482.xml')
ec_bu = read_sbml_model('../models/ec/ec_Bacteroides_uniformis_ATCC_8492.xml')
ec_ec = read_sbml_model('../models/ec/ec_Escherichia_coli_ED1a.xml')
ec_fn = read_sbml_model('../models/ec/ec_Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586.xml')
ec_ri = read_sbml_model('../models/ec/ec_Roseburia_intestinalis_L1_82.xml')
ec_sp = read_sbml_model('../models/ec/ec_Streptococcus_parasanguinis_ATCC_15912.xml')
ec_ss = read_sbml_model('../models/ec/ec_Streptococcus_salivarius_DSM_20560.xml')

# Abundances

abundances = {'M_Bacteroides_thetaiotaomicron_VPI_5482': 1,
              'M_Bacteroides_uniformis_ATCC_8492': 4,
              'M_Escherichia_coli_ED1a': 0.5,
              'M_Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586': 1,
              'M_Roseburia_intestinalis_L1_82': 1,
              'M_Streptococcus_parasanguinis_ATCC_15912': 1,
              'M_Streptococcus_salivarius_DSM_20560': 1}

# Media

M1 = {'R_EX_fe2[e]': (-0.000000033, 0)}
M2 = {'R_EX_fe2[e]': (0, 0), 'R_EX_pheme[e]': (0, 0)}
M3 = {'R_EX_fe2[e]': (-0.0000000082, 0),
      'R_EX_pheme[e]': (-0.0000000007, 0)
      }
M4 = {'R_EX_fe2[e]': (-0.000000033, 0),
      'R_EX_pheme[e]': (-0.0000000027, 0)
      }
M5 = {'R_EX_pheme[e]': (-0.0000000027, 0)}
M6 = {'R_EX_fe2[e]': (-0.000000125, 0)}
M7 = {'R_EX_pheme[e]': (-0.000000125, 0)}
M8 = {'R_EX_fe2[e]': (-0.000002, 0),
      'R_EX_pheme[e]': (-0.000000125, 0)
      }
M9 = {'R_EX_fe2[e]': (-0.000000125, 0),
      'R_EX_pheme[e]': (-0.000000125, 0)
      }
M10 = {}


def map_abundances(sample: List) -> List:
    '''
    Maps relative abundances to models in sample.

    '''
    sample_abun: List = []
    for model in sample:
        sample_abun.append(abundances[str(model.id)])
    return sample_abun


def adjust_media(media: dict, com: CommunityModel):
    model_ids_l = sorted(com.model_ids)
    for m in model_ids_l:
        media[f'R_prot_pool_exchange_{m}'] = (0, 1)
    return media


def generate_com(sample: List):
    # TODO: group samples by ec and non-ec
    com = CommunityModel(sample, abundances=map_abundances(sample))
    sim = com.get_community_model()
    sol_all: List = []
    sol_single: List = []
    model_ids_l = sorted(com.model_ids)
    model_ids: str = " + ".join(model_ids_l)
    media_l = [M1, M2, M3, M4, M5, M6, M7, M8, M9, M10]
    media_name_l = ['33 uM Fe', '0 uM Fe/H',
                    '0,7 uM H; 8,2 uM Fe', '2,7 uM H; 33 uM Fe',
                    '2,7 uM H', '125 uM Fe',
                    '125 uM H', '125 uM H; 2 mM Fe',
                    '125 uM H; 125 uM Fe', 'No constraints']
    for m in media_l:
        m = adjust_media(m, com)
        print(m)
        res: SimulationResult = sim.simulate(constraints=m)
        biomass_all = res.objective_value
        biomass_single = ', '.join(str(num) for num in res.find('R_biomass')['Flux rate'].to_list())
        sol_all.append(biomass_all)
        sol_single.append(biomass_single)
    data_all = {'Constraints': media_name_l, model_ids: sol_all}
    data_single = {'Constraints': media_name_l, model_ids: sol_single}
    df_all = pd.DataFrame(data_all)
    df_all = df_all.set_index('Constraints')
    df_single = pd.DataFrame(data_single)
    df_single = df_single.set_index('Constraints')
    return df_all, df_single


def main():
    merge_all: List = []
    merge_single: List = []
    with ProcessPoolExecutor() as p:
        for result_all, result_single in p.map(generate_com, [[ec_bt, ec_bu], [ec_bt, ec_ec], [ec_bt, ec_fn], [ec_bt, ec_ri], [ec_bt, ec_sp]]):
            merge_all.append(result_all)
            merge_single.append(result_single)
    df_all = pd.concat(merge_all, axis=1)
    print('All:\n')
    print(df_all)
    df_single = pd.concat(merge_single, axis=1)
    print('Single:\n')
    print(df_single)
    df_all.to_csv('../data/results/ec_fba_all.csv')
    df_single.to_csv('../data/results/ec_fba_single.csv')


if __name__ == "__main__":
    main()
