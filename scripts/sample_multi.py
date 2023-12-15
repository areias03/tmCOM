from cobra.io import read_sbml_model
from mewpy.cobra.com import *
from mewpy.model.com import CommunityModel
from mewpy.simulation import Environment 

import numpy as np
import pandas as pd

import sys
import itertools
from concurrent.futures import ProcessPoolExecutor

#Abundances

abundances = {'bt': 1, 'bu': 4, 'cc': 0.25, 'ec': 0.5,'fn': 1,'ri': 1,'sp': 1,'ss': 1}

#EC models

ec_bt = read_sbml_model('../models/ec/ec_Bacteroides_thetaiotaomicron_VPI_5482.xml')
ec_bu = read_sbml_model('../models/ec/ec_Bacteroides_uniformis_ATCC_8492.xml')
ec_ec = read_sbml_model('../models/ec/ec_Escherichia_coli_ED1a.xml')
ec_fn = read_sbml_model('../models/ec/ec_Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586.xml')
ec_ri = read_sbml_model('../models/ec/ec_Roseburia_intestinalis_L1_82.xml')
ec_sp = read_sbml_model('../models/ec/ec_Streptococcus_parasanguinis_ATCC_15912.xml')
ec_ss = read_sbml_model('../models/ec/ec_Streptococcus_salivarius_DSM_20560.xml')
ec_cc = read_sbml_model('../models/ec/ec_Coprococcus_comes_ATCC_27758.xml')

#Non-EC models

bt = read_sbml_model('../models/non-ec/agora/Bacteroides_thetaiotaomicron_VPI_5482.xml')
bu = read_sbml_model('../models/non-ec/agora/Bacteroides_uniformis_ATCC_8492.xml')
ec = read_sbml_model('../models/non-ec/agora/Escherichia_coli_ED1a.xml')
fn = read_sbml_model('../models/non-ec/agora/Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586.xml')
ri = read_sbml_model('../models/non-ec/agora/Roseburia_intestinalis_L1_82.xml')
sp = read_sbml_model('../models/non-ec/agora/Streptococcus_parasanguinis_ATCC_15912.xml')
ss = read_sbml_model('../models/non-ec/agora/Streptococcus_salivarius_DSM_20560.xml')
cc = read_sbml_model('../models/non-ec/agora/Coprococcus_comes_ATCC_27758.xml')


def main(arguments: Tuple[int, bool, str]) -> None:

    sample,enz,cons = arguments

    print(f'Conditions: Sample {sample}, Enzimatic constrained - {enz}, Conditions - {cons}')

    sample_name = 'sample' + str(sample)

     
    if sample == 1 and enz == True:
        community = CommunityModel([ec_bt,ec_bu,ec_ec],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ec']])
    elif sample == 1 and not enz:
        community = CommunityModel([bt,bu,ec],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ec']])
    elif sample == 2 and enz:
        community = CommunityModel([ec_bt,ec_bu,ec_fn],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['fn']])
    elif sample == 2 and not enz:
        community = CommunityModel([bt,bu,fn],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['fn']])
    elif sample == 3 and enz:
        community = CommunityModel([ec_bt,ec_bu,ec_ri],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ri']])
    elif sample == 3 and not enz:
        community = CommunityModel([bt,bu,ri],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ri']])
    elif sample == 4 and enz:
        community = CommunityModel([ec_bt,ec_bu,ec_sp],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['sp']])
    elif sample == 4 and not enz:
        community = CommunityModel([bt,bu,sp],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['sp']])
    elif sample == 5 and enz:
        community = CommunityModel([ec_bt,ec_bu,ec_ss],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ss']])
    elif sample == 5 and not enz:
        community = CommunityModel([bt,bu,ss],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ss']])
    elif sample == 6 and enz:
        community = CommunityModel([ec_bt,ec_bu,ec_cc], flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['cc']])
    elif sample == 6 and not enz:
        community = CommunityModel([bt,bu,cc], flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['cc']])
    else:
        sys.exit('No possible combination found!')

    if enz:
        sample_name = 'ec_' + sample_name

    if cons == 'Low Iron':
        sample_name = str(sample_name) + '_li'

    print(f'Writing results to ../data/results/{sample_name}.txt')

    with open(f'../data/results/{sample_name}.txt','w') as f:

        f.write(f'Conditions: Sample {sample}, Enzimatic constrained - {enz}, Conditions - {cons}\n')
        f.write('\n')

        f.write('##################################\n')
        f.write('COMMUNITY MODEL CREATION\n')
        f.write('##################################\n')

        model_ids = sorted(community.model_ids)

        f.write(f'Models in the sample:\n')
        for model in model_ids:
            f.write(f'{model}\n')
        f.write('\n')


        sim = community.get_community_model()

        f.write('\n')
        f.write(f'Number of reactions in the community model: {len(sim.reactions)}\n')
        f.write('\n')


        M9 = Environment.from_model(bu)

        f.write(f'Environment from model: {bu.id}\n')
        f.write('\n')
        print(M9,file=f)
        f.write('\n')

        f.write('##################################\n')
        f.write('SIMULATION\n')
        f.write('##################################\n')

        sim.set_environmental_conditions(M9)
       

        if enz and cons =='Low Iron':
            constraints = {'R_EX_fe2[e]':(-0.00007,0),'R_prot_pool_exchange_M_Bacteroides_thetaiotaomicron_VPI_5482':(0,1),'R_prot_pool_exchange_M_Bacteroides_uniformis_ATCC_8492':(0,1),f'R_prot_pool_exchange_{model_ids[2]}':(0,1)}
        elif enz and cons =='Default':
            constraints = {'R_prot_pool_exchange_M_Bacteroides_thetaiotaomicron_VPI_5482':(0,1),'R_prot_pool_exchange_M_Bacteroides_uniformis_ATCC_8492':(0,1),f'R_prot_pool_exchange_{model_ids[2]}':(0,1)}
        elif not enz and cons =='Low Iron':
            constraints = {'R_EX_fe2[e]':(-0.00007,0)}
        else:
            constraints = {}
        
        f.write('Constraints\n')
        f.write(pd.DataFrame.from_dict(constraints).to_string())
        f.write('\n')

        solution = sim.simulate(constraints=constraints)

        f.write('FBA - Flux Balance Analysis\n')

        f.write(f'Community growth: {solution.objective_value}\n')

        f.write('Individual growth\n')
        f.write(solution.find('R_biomass', sort=True,show_nulls=True).to_string())
        f.write('\n')

        f.write('Exchange reactions and fluxes\n')
        f.write(solution.find('R_EX').to_string())
        f.write('\n')

        f.write('Iron exchange reactions and fluxes\n')
        f.write(solution.find('R_EX_fe2').to_string())
        f.write('\n')

        f.write('SteadyCom\n')
        try:
            solution = SteadyCom(community, constraints=constraints)
            print(solution,file= f)
        except ZeroDivisionError:
            f.write('No possible solution.\n')
        f.write('\n')

        f.write('Cross-feeding interactions:\n')
        f.write(solution.cross_feeding(as_df=True).dropna().sort_values('rate', ascending=False).to_string())
        f.write('\n')

        

        f.write('SteadyCom - Variability Analysis\n')
        f.write('\n')
              
        l_va = np.linspace(0.1,1.0,num=10)

        for va in l_va:
            va = round(va, 1)
            variability = SteadyComVA(community, obj_frac=va, constraints=M9)
            f.write(f'Strain\tMin\tMax\tVariability - {va}\n')
            for strain, (lower, upper) in variability.items():
                f.write(f'{strain}\t{lower:.1%}\t{upper:.1%}\n')
       

        f.write('\n')
        f.write('##################################\n')
        f.write('COMMUNITY ANALYSIS\n')
        f.write('##################################\n')

        f.write('SMETANA - Species Metabolic Interaction Analysis\n')

        f.write('SCS (species coupling score):\n')
        SCS = sc_score(community)
        f.write(pd.DataFrame.from_dict(SCS).to_string())
        f.write('\n')
              
        f.write('MUS (metabolite uptake score):\n')
        MUS = mu_score(community)
        f.write(pd.DataFrame.from_dict(MUS).to_string())
        f.write('\n')

        
        f.write('MPS (metabolite production score):\n')
        MPS = mp_score(community,environment=M9)
        f.write(pd.DataFrame.from_dict(MPS).to_string())
        f.write('\n')
        
              
        f.write('MRO (metabolic resource overlap):\n')
        score, MRO = mro_score(community,environment=M9)
        f.write(f'Community score: {score}\n')

        f.write('Total competition for resources:\n')
        f.write(pd.DataFrame.from_dict(MRO.community_medium).to_string())
        f.write('')
        f.write('By individual:\n')

        for ind in MRO.individual_media.keys():
            f.write(f'Strain:{ind}\t{", ".join(met for met in MRO.individual_media[ind])}\n\n')
        f.close()
    


if __name__ == "__main__":

    sample: List[int] = list(range(1, 7))
    enz: List[bool] = [True, False]
    li: List[str] = ["Default", "Low Iron"]

    args: List[Tuple[int, bool, str]] = list(itertools.product(sample, enz, li))

    print('Possible combinations:', f"Sample: {sample} | Enz. constraints: {enz} | Conditions: {li}")

    with ProcessPoolExecutor() as pool:
        pool.map(main,args)
