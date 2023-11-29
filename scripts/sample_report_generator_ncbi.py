from cobra.io import read_sbml_model
from mewpy import get_simulator
from mewpy.cobra.com import *
from com import CommunityModel
from mewpy.simulation import Environment
from mewpy.cobra.com import SteadyComVA

import numpy as np
from multiprocessing.pool import ThreadPool
import time
st = time.time()

bt = read_sbml_model('../models/non-ec/ncbi_refseq/vpi5482.xml')
simbt = get_simulator(bt)
simbt.set_objective("Growth")

bu = read_sbml_model('../models/non-ec/ncbi_refseq/atcc8492.xml')
simbu = get_simulator(bu)
simbu.set_objective("Growth")

ec = read_sbml_model('../models/non-ec/ncbi_refseq/ed1a.xml')
simec = get_simulator(ec)
simec.set_objective("Growth")
print()
print()

print('##################################')
print('COMMUNITY MODEL CREATION')
print('##################################\n')

community = CommunityModel([bt,bu,ec],flavor='cobra')


sim = community.get_community_model()

print(f'Number of reactions in the community model: {len(sim.reactions)}\n')

from mewpy.simulation import Environment
M9 = Environment.from_model(ec)

print('Environment from model\n')
print(M9)
print()

print('##################################')
print('SIMULATION')
print('##################################\n')

solution = sim.simulate(constraints=M9)

print('FBA - Flux Balance Analysis\n')

print(f'Community growth: {solution.objective_value}\n')

print('Individual growth\n')
print(solution.find('biomass', sort=True,show_nulls=True))
print()

print('Exchange reactions and fluxes\n')
print(solution.find('EX'))
print()

print('SteadyCom\n')

solution = SteadyCom(community, constraints=M9)
print(solution)
print()

print('Cross-feeding interactions:\n')
print(solution.cross_feeding(as_df=True).dropna().sort_values('rate', ascending=False))
print()
      
print('SteadyCom - Variability Analysis\n')
      
l_va = np.linspace(0.1,1.0,num=10)

for va in l_va:
    va = round(va, 1)
    variability = SteadyComVA(community, obj_frac=va, constraints=M9)
    print(f'Strain\tMin\tMax\tVariability - {va}')
    for strain, (lower, upper) in variability.items():
        print(f'{strain}\t{lower:.1%}\t{upper:.1%}')

print()
print('##################################')
print('COMMUNITY ANALYSIS')
print('##################################\n')
      
print('Jacard Similarity Matrices\n')

mets, rxns, over = jaccard_similarity_matrices([glc_ko, nh4_ko,wildtype])

print('Metabolite overlap\n')
print(mets)
print()
mets_html = mets.to_html("../data/images/metsover_sample_test.html")

print('Reactions overlap\n')
print(rxns)
print()
rxns_html = rxns.to_html("../data/images/rxnsover_sample_test.html")

print('Uptake overlap\n')
print(over)
print()
over_html = over.to_html("../data/images/uptakeover_sample_test.html")
print('SMETANA - Species Metabolic Interaction Analysis\n')

print('SCS (species coupling score):\n')
print(sc_score(community))
print()
      
print('MUS (metabolite uptake score):\n')
MUS = mu_score(community)
print(pd.DataFrame.from_dict(MUS))
print()

print('MPS (metabolite production score):\n')
MPS = mp_score(community,environment=M9)
print(pd.DataFrame.from_dict(MPS))
print()
      
print('MRO (metabolic resource overlap):\n')
score, MRO = mro_score(community,environment=M9)
print(f'Community score: {score}\n')

print('Total competition for resources:\n')
print(MRO.community_medium)
print()
print('By individual:\n')

for ind in MRO.individual_media.keys():
    print(f'Strain:{ind}\t{", ".join(met for met in MRO.individual_media[ind])}\n\n')
    
et = time.time()

elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')