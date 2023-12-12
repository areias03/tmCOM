from cobra.io import read_sbml_model
from mewpy.cobra.com import *
from mewpy.model.com import CommunityModel
from mewpy.simulation import Environment 

import numpy as np
import pandas as pd

from concurrent.futures import ThreadPoolExecutor

#Non-EC models

bt = read_sbml_model('../models/non-ec/agora/Bacteroides_thetaiotaomicron_VPI_5482.xml')

bu = read_sbml_model('../models/non-ec/agora/Bacteroides_uniformis_ATCC_8492.xml')

ec = read_sbml_model('../models/non-ec/agora/Escherichia_coli_ED1a.xml')

fn = read_sbml_model('../models/non-ec/agora/Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586.xml')

ri = read_sbml_model('../models/non-ec/agora/Roseburia_intestinalis_L1_82.xml')

sp = read_sbml_model('../models/non-ec/agora/Streptococcus_parasanguinis_ATCC_15912.xml')

ss = read_sbml_model('../models/non-ec/agora/Streptococcus_salivarius_DSM_20560.xml')

#EC models

ec_bt = read_sbml_model('../models/ec/ec_Bacteroides_thetaiotaomicron_VPI_5482.xml')

ec_bu = read_sbml_model('../models/ec/ec_Bacteroides_uniformis_ATCC_8492.xml')

ec_ec = read_sbml_model('../models/ec/ec_Escherichia_coli_ED1a.xml')

ec_fn = read_sbml_model('../models/ec/ec_Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586.xml')

ec_ri = read_sbml_model('../models/ec/ec_Roseburia_intestinalis_L1_82.xml')

ec_sp = read_sbml_model('../models/ec/ec_Streptococcus_parasanguinis_ATCC_15912.xml')

ec_ss = read_sbml_model('../models/ec/ec_Streptococcus_salivarius_DSM_20560.xml')

#Abundances

abundances = {'bt': 1, 'bu': 4, 'ec': 0.5,'fn': 1,'ri': 1,'sp': 1,'ss': 1}

#Non-ec samples

sample1 = CommunityModel([bt,bu,ec],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ec']])
sample2 = CommunityModel([bt,bu,fn],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['fn']])
sample3 = CommunityModel([bt,bu,ri],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ri']])
sample4 = CommunityModel([bt,bu,sp],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['sp']])
sample5 = CommunityModel([bt,bu,ss],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ss']])

#EC samples

ec_sample1 = CommunityModel([ec_bt,ec_bu,ec_ec],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ec']])
ec_sample2 = CommunityModel([ec_bt,ec_bu,ec_fn],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['fn']])
ec_sample3 = CommunityModel([ec_bt,ec_bu,ec_ri],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ri']])
ec_sample4 = CommunityModel([ec_bt,ec_bu,ec_sp],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['sp']])
ec_sample5 = CommunityModel([ec_bt,ec_bu,ec_ss],flavor='reframed', add_compartments=True, merge_biomasses=True,abundances=[abundances['bt'],abundances['bu'],abundances['ss']])


def main(sample:int,cons:str = 'Default',ec:bool = False):

    if sample == 1 and ec:
        community = ec_sample1
        sample_name = 'ec_sample1'
    elif sample == 1 and not ec:
        community = sample1
        sample_name = 'sample1'
    elif sample == 2 and ec:
        community = ec_sample2
        sample_name = 'ec_sample2'
    elif sample == 2 and not ec:
        community = sample2
        sample_name = 'sample2'
    elif sample == 3 and ec:
        community = ec_sample3
        sample_name = 'ec_sample3'
    elif sample == 3 and not ec:
        community = sample3
        sample_name = 'sample3'
    elif sample == 4 and ec:
        community = ec_sample4
        sample_name= 'ec_sample4'
    elif sample == 4 and not ec:
        community = sample4
        sample_name= 'sample4'
    elif sample == 5 and ec:
        community = ec_sample5
        sample_name= 'ec_sample5'
    else:
        community = sample5
        sample_name= sample5

    with open(f'../data/results/{sample_name}.txt','w') as f:

        f.write('##################################\n')
        f.write('COMMUNITY MODEL CREATION\n')
        f.write('##################################\n')


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
       
        model_ids = sorted(community.model_ids)

        if ec and cons =='Low iron':
            constraints = {'R_EX_fe2(e)':(-0.00007,0),'R_prot_pool_exchange_M_Bacteroides_thetaiotaomicron_VPI_5482':(0,1),'R_prot_pool_exchange_M_Bacteroides_uniformis_ATCC_8492':(0,1),f'R_prot_pool_exchange_{model_ids[2]}':(0,1)}
        elif ec and cons =='Default':
            constraints = {'R_prot_pool_exchange_M_Bacteroides_thetaiotaomicron_VPI_5482':(0,1),'R_prot_pool_exchange_M_Bacteroides_uniformis_ATCC_8492':(0,1),f'R_prot_pool_exchange_{model_ids[2]}':(0,1)}
        elif not ec and cons =='Low iron':
            constraints = {'R_EX_fe2(e)':(-0.00007,0)}
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

        '''

        f.write('SteadyCom - Variability Analysis\n')
              
        l_va = np.linspace(0.1,1.0,num=10)

        for va in l_va:
            va = round(va, 1)
            variability = SteadyComVA(community, obj_frac=va, constraints=M9)
            f.write(f'Strain\tMin\tMax\tVariability - {va}')
            for strain, (lower, upper) in variability.items():
                f.write(f'{strain}\t{lower:.1%}\t{upper:.1%}')
        '''
        f.write('\n')
        f.write('##################################\n')
        f.write('COMMUNITY ANALYSIS\n')
        f.write('##################################\n')
              
        f.write('Jacard Similarity Matrices\n')

        mets, rxns, over = jaccard_similarity_matrices([bt,bu,ec])

        f.write('Metabolite overlap\n')
        f.write(mets.to_string())
        f.write('\n')

        f.write('Reactions overlap\n')
        f.write(rxns.to_string())
        f.write('\n')

        f.write('Uptake overlap\n')
        f.write(over.to_string())
        f.write('\n')

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
    with ThreadPoolExecutor() as pool:
        pool.map(main,zip((1,2,3,4,5),(True,False),('Default','Low iron')))
