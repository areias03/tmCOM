from cobra.io import read_sbml_model
from mewpy.cobra.com import Simulator, SteadyCom, SteadyComVA, mu_score, mp_score, mip_score, mro_score, sc_score, minimal_medium, build_problem
from mewpy.model.com import CommunityModel
from mewpy.solvers.solution import SimulationResult

import numpy as np
import pandas as pd

from typing import List
from concurrent.futures import ProcessPoolExecutor


# EC models

ec_bt = read_sbml_model('../../models/ec/ec_Bacteroides_thetaiotaomicron_VPI_5482.xml')
ec_bu = read_sbml_model('../../models/ec/ec_Bacteroides_uniformis_ATCC_8492.xml')
ec_fn = read_sbml_model('../../models/ec/ec_Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586.xml')
ec_ri = read_sbml_model('../../models/ec/ec_Roseburia_intestinalis_L1_82.xml')
ec_sp = read_sbml_model('../../models/ec/ec_Streptococcus_parasanguinis_ATCC_15912.xml')


prot_pool = {'M_Bacteroides_thetaiotaomicron_VPI_5482': (0,0.075),
              'M_Bacteroides_uniformis_ATCC_8492': (0,0.075),
              'M_Escherichia_coli_ED1a': (0,0.13),
              'M_Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586': (0,0.046),
              'M_Roseburia_intestinalis_L1_82': (0,0.078),
              'M_Streptococcus_parasanguinis_ATCC_15912': (0,0.082),
              'M_Streptococcus_salivarius_DSM_20560': (0,0.048)}


# Media

constraints = {'R_EX_glc_D[e]':(-5,100), # D-glucose
               'R_EX_fru[e]':(-1,0), # Fructose
               'R_EX_cellb[e]':(-1,100), # Cellobiose
               'R_EX_malt[e]':(-1,100), # Maltose
               'R_EX_lcts[e]':(-1,100), # Lactose
               'R_EX_inulin[e]':(-1,100), # Inulin
               'R_EX_nad[e]':(-1,100), # NAD
               'R_EX_his_L[e]':(-0.171,100), # L-Histidine
               'R_EX_ile_L[e]':(-0.24,100), # L-Isoleucine
               'R_EX_leu_L[e]':(-1,100), # L-Leucine
               'R_EX_met_L[e]':(-0.125,100), # L-Methionine
               'R_EX_val_L[e]':(-0.7,100), # L-Valine
               'R_EX_arg_L[e]':(-0.72,100), # L-Arginine
               'R_EX_cys_L[e]':(-0.5,100), # L-Cysteine
               'R_EX_glu_L[e]':(-0.6,100), # L-Glutamate
               'R_EX_phe_L[e]':(-0.4,100), # L-Phenylalanine
               'R_EX_pro_L[e]':(-0.7,100), # L-Proline
               'R_EX_asn_L[e]':(-0.5,100), # L-Asparagine
               'R_EX_asp_L[e]':(-0.05,100), # L-Aspartate
               'R_EX_gln_L[e]':(-0.6,100), # L-Glutamine
               'R_EX_ser_L[e]':(-0.5,100), # L-Serine
               'R_EX_thr_L[e]':(-0.5,100), # L-Threonine
               'R_EX_ala_L[e]':(-0.4,100), # L-Alanine
               'R_EX_gly[e]':(-0.3,100), # Glycine
               'R_EX_lys_L[e]':(-0.5,100), # L-Lysine
               'R_EX_trp_L[e]':(-0.2,100), # L-Tryptophan
               'R_EX_tyr_L[e]':(-0.3,100), # L-Tyrosine
               'R_EX_ade[e]':(-0.011,100), # Adenine
               'R_EX_gua[e]':(-0.0056,100), # Guanine
               'R_EX_ura[e]':(-0.023,100), # Uracil
               'R_EX_xan[e]':(-0.0038,100), # Xanthine
               'R_EX_mgcl2[e]':(-0.386,100), # Magnesium Chloride
               'R_EX_btn[e]':(-6.02,100), # Biotin
               'R_EX_ribflv[e]':(-0.9,100), # Riboflavin
               'R_EX_fol[e]':(-0.56,100), # Folate
               'R_EX_nac[e]':(-0.9,100), # Niacin
               'R_EX_inost[e]':(-0.2,100), # Inositol
               'R_EX_gthrd[e]':(-0.9,100), # L-Glutathione reduced
}


M1 = {'R_EX_fe2[e]': (-0.000000033, 0), 'R_EX_pheme[e]': (0, 0)}
M2 = {'R_EX_fe2[e]': (0, 0), 'R_EX_pheme[e]': (0, 0)}
M3 = {'R_EX_fe2[e]': (-0.0000000082, 0),
      'R_EX_pheme[e]': (-0.0000000007, 0)
      }
M4 = {'R_EX_fe2[e]': (-0.000000033, 0),
      'R_EX_pheme[e]': (-0.0000000027, 0)
      }
M5 = {'R_EX_pheme[e]': (-0.0000000027, 0), 'R_EX_fe2[e]': (0, 0)}
M6 = {'R_EX_fe2[e]': (-0.000000125, 0), 'R_EX_pheme[e]': (0, 0)}
M7 = {'R_EX_pheme[e]': (-0.000000125, 0), 'R_EX_fe2[e]': (0, 0)}
M8 = {'R_EX_fe2[e]': (-0.000002, 0),
      'R_EX_pheme[e]': (-0.000000125, 0)
      }
M9 = {'R_EX_fe2[e]': (-0.000000125, 0),
      'R_EX_pheme[e]': (-0.000000125, 0)
      }
M10 = {}


def adjust_media(media: dict, com: CommunityModel):
    model_ids_l = sorted(com.model_ids)
    for m in model_ids_l:
        media[f'R_prot_pool_exchange_{m}'] = prot_pool[m]
    return media


def generate_com(sample: List):
    com = CommunityModel(sample)
    sim = com.get_community_model()
    solver = build_problem(com)
    sol_all: List = []
    sol_single: List = []
    sol_sc: List = []
    sol_abun_sc: List = []
    model_ids_l = sorted(com.model_ids)
    model_ids: str = " + ".join(model_ids_l)
    print(model_ids)
    media_l = [M1, M2, M3, M4, M5, M6, M7, M8, M9, M10]
    media_name_l = ['33 uM Fe', '0 uM Fe/H',
                    '0,7 uM H; 8,2 uM Fe', '2,7 uM H; 33 uM Fe',
                    '2,7 uM H', '125 uM Fe',
                    '125 uM H', '125 uM H; 2 mM Fe',
                    '125 uM H; 125 uM Fe', 'No constraints']
    for m_ind, m in enumerate(media_l):
        m = m | constraints
        m_pp = m.copy()
        m_pp = adjust_media(m_pp, com)
        res: SimulationResult = sim.simulate(constraints=m_pp)
        biomass_all = res.objective_value
        try:
            biomass_bt = res.find('R_biomass_M_Bacteroides_thetaiotaomicron_VPI_5482')['Flux rate'].to_list()[0]
        except IndexError:
            biomass_bt = 0
        try:
            res_sc = SteadyCom(com,constraints=m,solver=solver)
            growth_sc = res_sc.growth
            abun_sc = list(res_sc.abundance.values())[0]
        except ZeroDivisionError:
            growth_sc = 0
            abun_sc = 0           
        data_single = [media_name_l[m_ind], model_ids_l[1], biomass_bt]
        data_abun_sc = [media_name_l[m_ind], model_ids_l[1], abun_sc]
        sol_all.append(biomass_all)
        sol_sc.append(growth_sc)
        sol_single.append(data_single)
        sol_abun_sc.append(data_abun_sc)
    data_all = {'Constraints': media_name_l, model_ids: sol_all}
    data_sc = {'Constraints': media_name_l, model_ids: sol_sc}
    df_all = pd.DataFrame(data_all)
    df_all = df_all.set_index('Constraints')
    df_sc = pd.DataFrame(data_sc)
    df_sc = df_sc.set_index('Constraints')
    df_single = pd.DataFrame(sol_single)
    df_single.columns = ['Constraints','Partner','Growth']
    df_abun_sc = pd.DataFrame(sol_abun_sc)
    df_abun_sc.columns = ['Constraints','Partner','Growth']
    return df_all, df_single, df_sc, df_abun_sc


def main():
    merge_all: List = []
    merge_single: List = []
    merge_sc: List = []
    merge_abun_sc: List = []
    with ProcessPoolExecutor() as p:
        for result_all, result_single, result_sc, result_abun_sc in p.map(generate_com, [[ec_bt,ec_bu], [ec_bt,ec_fn], [ec_bt, ec_ri], [ec_bt, ec_sp]]):
            merge_all.append(result_all)
            merge_single.append(result_single)
            merge_sc.append(result_sc)
            merge_abun_sc.append(result_abun_sc)
    df_all = pd.concat(merge_all, axis=1)
    print('FBA:\n\n')
    print('All:\n')
    print(df_all)
    df_single = pd.concat(merge_single, axis=0)
    print('Single:\n')
    print(df_single)
    df_sc = pd.concat(merge_sc, axis=1)
    print('SteadyCom:\n\n')
    print('All:\n')
    print(df_sc)
    df_abun_sc = pd.concat(merge_abun_sc, axis=0)
    print('Single:\n')
    print(df_abun_sc)
    df_all.to_csv('results/ec_fba_all.csv')
    df_single.to_csv('results/ec_fba_single_all.csv')
    df_sc.to_csv('results/ec_sc_all.csv')
    df_abun_sc.to_csv('results/ec_sc_single_all.csv')


if __name__ == "__main__":
    main()
