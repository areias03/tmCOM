from tqdm.auto import tqdm
from reframed.io.sbml import load_cbmodel,save_cbmodel
from cobra.io import read_sbml_model, write_sbml_model
from mewpy.cobra.util import add_enzyme_constraints,convert_gpr_to_dnf,split_isozymes
from mewpy.simulation import get_simulator
from mewpy.simulation.environment import Environment as Environment
from mewpy.util.request import retreive_gene,retreive_protein,get_smiles,brenda_query
import pandas as pd
import numpy as np
from urllib.request import urlopen
from functools import reduce
import json
from ast import literal_eval
import requests

logger.info("Loading SBML model and scraping annotations...")

filepath = "../models/non-ec/Bacteroides_uniformis_ATCC_8492.xml"
model = read_sbml_model(filepath)
model2 = load_cbmodel(filepath)
organism = 'Bacteroides thetaiotaomicron'

sim = get_simulator(model)
sim.set_objective("biomass")

sim2 = get_simulator(model2)
sim2.set_objective("R_biomass")

ls_ge = []

for ge in sim.genes:
    i = sim.genes.index(ge)
    gene = sim2.genes[i]
    rxns = sim.get_gene(ge).reactions
    for rx in rxns:
        anno = sim.get_reaction(rx)['annotations']
        seed_id = anno.get('seed.reactions')
        rxn_name = sim.get_reaction(rx).name
        ecnumber = anno.get('ec-code')
        metanetx = anno.get('metanetx.reaction')
        kegg = anno.get('kegg.reaction')
        res = [gene,rx,rxn_name,seed_id,metanetx,kegg,ecnumber]
        ls_ge.append(res)

df_ge = pd.DataFrame(ls_ge,columns=[['Gene','Reaction','Name','ModelSEED_id','MetaNetX','KEGG_id','ecNumber']])

seed_id = df_ge['ModelSEED_id'].values.tolist()

seed_id = [reduce(lambda x: x, inner_list) for inner_list in seed_id]

metanetx_id = df_ge['MetaNetX'].values.tolist()

metanetx_id = [reduce(lambda x: x, inner_list) for inner_list in metanetx_id]

kegg_id = df_ge['KEGG_id'].values.tolist()

kegg_id = [reduce(lambda x: x, inner_list) for inner_list in kegg_id]

logger.info("Querying MODELSEED database...")

SOLR_URL='https://modelseed.org'
ls_name = []
ls_kegg = []
ls_bigg = []

for mseed_id in tqdm(seed_id):
    i = seed_id.index(mseed_id)
    if mseed_id == None:
        ls_name.append('None')
        ls_kegg.append('None')
        ls_bigg.append('None')
    else:
        try:
            connection = urlopen(SOLR_URL+f'/solr/reactions/select?wt=json&q=id:{mseed_id}&fl=name,id,formula,charge,aliases')
            response = json.load(connection)
            for document in response['response']['docs']:  
                ms_name = document.get('name')
                ls_alias = document.get('aliases')
                ms_bigg = list(filter(lambda a: 'BiGG:' in a, document.get('aliases')))
                ms_kegg = list(filter(lambda a: 'KEGG:' in a, document.get('aliases')))
                if len(ms_bigg)== 0 and len(ms_kegg)== 0:
                    ms_bigg = 'None'
                    ms_kegg = 'None'
                elif len(ms_bigg)== 0 and len(ms_kegg)!= 0:
                    ms_bigg = 'None'
                    ms_kegg = list(ms_kegg)[0]
                    ms_kegg = ms_kegg.replace('KEGG: ','')
                elif len(ms_bigg)!= 0 and len(ms_kegg)== 0:
                    ms_kegg = 'None'
                    ms_bigg = list(ms_bigg)[0]
                    ms_bigg = ms_bigg.replace('BiGG: ','')
                else:
                    ms_kegg = list(ms_kegg)[0]
                    ms_kegg = ms_kegg.replace('KEGG: ','')
                    ms_bigg = list(ms_bigg)[0]
                    ms_bigg = ms_bigg.replace('BiGG: ','')    
                ls_name.append(ms_name)
                ls_bigg.append(ms_bigg)
                ls_kegg.append(ms_kegg)           
        except:
            ls_name.append('None')
            ls_kegg.append('None')
            ls_bigg.append('None')
            
new_kegg = [next(filter(None, i)) for i in zip(ls_kegg, kegg_id)]

df_ge['BIGG_id'] = ls_bigg
df_ge['KEGG_id'] = new_kegg

logger.info("Querying BIGG database...")

ls_bigg = df_ge['BIGG_id'].values.tolist()
bigg_ls = []

for bigg in tqdm(ls_bigg):
    for bigg_n in bigg:
        bigg_n = str(bigg_n)
        bigg_n = bigg_n.split(';')
        bigg_n = [x.strip(' ') for x in bigg_n]
        sub_bigg_ls = []
        for bi in bigg_n:
            if bi == 'None':
                pass
            else:
                url =f'http://bigg.ucsd.edu/api/v2/universal/reactions/{bi}'
                with requests.request("GET", url) as resp:
                    try:
                        resp.raise_for_status()  # raises exception when not a 2xx response
                        if resp.status_code != 204:
                            data = dict(resp.json())
                            ec_l = data['database_links']
                            if ec_l == None:
                                sub_bigg_ls.append(None)
                            else:
                                ec = [i['id'] for i in ec_l['EC Number']]
                                if ec == None:
                                    sub_bigg_ls.append(None)
                                sub_bigg_ls.append(ec)
                        else: 
                            sub_bigg_ls.append(None)
                    except:
                        sub_bigg_ls.append(None)
    bigg_ls.append(sub_bigg_ls)
    
ec_l = df_ge['ecNumber'].values.tolist()
new_l = [next(filter(None, i)) for i in zip(bigg_ls, ec_l)]
df_ge['ecNumber'] = new_l

logger.info("Obtaining reaction substrates and substrate SMILES...")

rx_l = df_ge['Reaction'].values.tolist()
ls_sub = []

for rxn in rx_l:
    for rx in rxn:
        sub = list(sim.get_substrates(rx).keys())
        ls_sub.append(sub)
    
df_ge["Substrates"] = ls_sub   

sub_na = df_ge['Substrates'].values.tolist()

ls_sub = []
ls_smile =[]



for sub_l in tqdm(sub_na):
    sub_ls_sub = []
    sub_ls_smile = []
    for sub_s in sub_l:
        for sub in sub_s:
            sub_name = sim.get_metabolite(sub).get('name')
            if "2,3-C" in sub_name:
                sub_name = sub_name.replace("2,3","2',3'")
            elif "2,3-c" in sub_name:
                sub_name = sub_name.replace("2,3","2',3'")
            elif "3-triphosphate" in sub_name:
                sub_name = sub_name.replace("3","3'")
            else:
                pass
            smile = get_smiles(sub_name)
            sub_ls_smile.append(smile)
            sub_ls_sub.append(sub_name)
        ls_sub.append(sub_ls_sub)
        ls_smile.append(sub_ls_smile)

df_ge['Substrate Name'] = ls_sub
df_ge['Substrate SMILES'] = ls_smile


ec_l = df_ge['ecNumber'].values.tolist()
ec_nl = []
ei = []

for es in ec_l:
    es = str(es)
    es = es.split(',')
    ei = []
    for sublist in es:
        sublist = str(sublist)
        sublist = sublist.strip("[[")
        sublist = sublist.strip(" ")
        sublist = sublist.strip("'")
        sublist = sublist.strip("[")
        sublist = sublist.strip("[")
        sublist = sublist.strip("'")
        sublist = sublist.strip("'")
        sublist = sublist.strip("]")
        sublist = sublist.strip("]]")  
        sublist = sublist.strip("'")
        if sublist not in ei:
            if sublist == 'None' and len(ei) > 0:
                pass
            elif '-' in sublist:
                pass
            else:
                ei.append(sublist)
    ec_nl.append(ei)

df_ge['ecNumber'] = ec_nl

logger.info("Querying BRENDA database...")

logger.info("Obtaining Kcats...")


from brendapyrser import BRENDA
from brendapyrser import EnzymePropertyDict


dataFile = '../../brenda_2023_1.txt'
brenda = BRENDA(dataFile)

kcat_ls = []
ec_ls = df_ge['ecNumber'].values.tolist()

for ec in tqdm(ec_ls):
    ec=str(ec)
    ec = ec.split(',')
    sub_kcat_ls = []
    #print(ec)
    for ec_n in ec:
        ec_n = ec_n.strip(',')
        ec_n = ec_n.strip('[')
        ec_n = ec_n.strip('[[')
        ec_n = ec_n.strip(']')
        ec_n = ec_n.strip(']]')
        ec_n = ec_n.strip(' ')
        ec_n = ec_n.strip('[')
        ec_n = ec_n.strip("'")
        ec_n = ec_n.strip(',')
        ec_n = ec_n.strip(' ')
        if "-" in ec_n:
            pass
        elif ec_n == None:
            sub_kcat_ls.append(None)
            continue
        else:
            try:
                r = brenda.reactions.get_by_id(ec_n)
                kcat_va = r.Kcatvalues.get_values()
                avg_kcat = sum(kcat_va)/len(kcat_va)
                sub_kcat_ls.append(avg_kcat)
            except:
                sub_kcat_ls.append(None)
    kcat_ls.append(sub_kcat_ls)

df_ge['Avg Kcat (by ec)'] = kcat_ls

kcat_ls = []
ec_ls = df_ge['ecNumber'].values.tolist()


for ec in tqdm(ec_ls):
    ec=str(ec)
    ec = ec.split(',')
    sub_kcat_ls = []
    for ec_n in ec:
        ec_n = ec_n.strip(',')
        ec_n = ec_n.strip('[')
        ec_n = ec_n.strip('[[')
        ec_n = ec_n.strip(']')
        ec_n = ec_n.strip(']]')
        ec_n = ec_n.strip(' ')
        ec_n = ec_n.strip('[')
        ec_n = ec_n.strip("'")
        if "-" in ec_n:
            pass
        else:
            try:
                r = brenda.reactions.get_by_id(ec_n)
                kcat_va = r.Kcatvalues.filter_by_organism(organism).get_values()
                avg_kcat = sum(kcat_va)/len(kcat_va)
                sub_kcat_ls.append(avg_kcat)
            except:
                sub_kcat_ls.append(None)
    kcat_ls.append(sub_kcat_ls)

df_ge['Avg Kcat (by ec and species)'] = kcat_ls

seq_ls = []
ec_ls = df_ge['ecNumber'].values.tolist()

logger.info("Obtaining protein sequences...")


for ec in tqdm(ec_ls):
    ec=str(ec)
    ec = ec.split(',')
    sub_seq_ls = []
    for ec_n in ec:
        ec_n = ec_n.strip(',')
        ec_n = ec_n.strip('[')
        ec_n = ec_n.strip('[[')
        ec_n = ec_n.strip(']')
        ec_n = ec_n.strip(']]')
        ec_n = ec_n.strip(' ')
        ec_n = ec_n.strip('[')
        ec_n = ec_n.strip("'")
        ec_n = ec_n.strip(' ')
        ec_n = ec_n.strip('[')
        if "-" in ec_n:
            pass
        else:
            try:
                from zeep import Client
                import hashlib

                wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
                password = hashlib.sha256("07042000Alex!".encode("utf-8")).hexdigest()
                client = Client(wsdl)
                parameters = ( "alexandreareias1718@gmail.com",password,f"ecNumber*{ec_n}","sequence*", "noOfAminoAcids*", "firstAccessionCode*", "source*", "id*", "organism*Bacteroides sp")
                resultString = client.service.getSequence(*parameters) 
                sub_seq_ls.append(resultString[0]['sequence'])
            except:
                sub_seq_ls.append(None)                  
    seq_ls.append(sub_seq_ls)
    
df_ge['Protein Sequence'] = seq_ls

logger.info("Obtaining molecular weights...")

mw_ls = []
ec_ls = df_ge['ecNumber'].values.tolist()


for ec in tqdm(ec_ls):
    ec=str(ec)
    ec = ec.split(',')
    sub_mw_ls = []
    for ec_n in ec:
        ec_n = ec_n.strip(',')
        ec_n = ec_n.strip('[')
        ec_n = ec_n.strip('[[')
        ec_n = ec_n.strip(']')
        ec_n = ec_n.strip(']]')
        ec_n = ec_n.strip(' ')
        ec_n = ec_n.strip('[')
        ec_n = ec_n.strip("'")
        ec_n = ec_n.strip(' ')
        ec_n = ec_n.strip('[')
        if "-" in ec_n:
            pass
        else:
            try:
                from zeep import Client
                import hashlib
                
                res_mw = 0
                wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
                password = hashlib.sha256("07042000Alex!".encode("utf-8")).hexdigest()
                client = Client(wsdl)
                parameters = ("alexandreareias1718@gmail.com",password,f"ecNumber*{ec_n}","molecularWeight*","molecularWeightMaximum*","commentary*","organism*","literature*" )
                resultString = client.service.getMolecularWeight(*parameters)
                for i in range(len(resultString)):
                    res_mw = res_mw + int(resultString[i]['molecularWeight'])
                res_mw = res_mw/len(resultString)
                sub_mw_ls.append(res_mw)
            except:
                sub_mw_ls.append(None)                  
    mw_ls.append(sub_mw_ls)

df_ge['Molecular Weight'] = mw_ls

logger.info("Preparing for DLKcat input...")

dk_prep = df_ge.drop(columns=['Substrates','ModelSEED_id','MetaNetX','KEGG_id','BIGG_id','ecNumber'])

dk_prep['Substrate Name'] = dk_prep['Substrate Name'].apply(literal_eval) #convert to list type
dk_prep['Substrate SMILES'] = dk_prep['Substrate SMILES'].apply(literal_eval) #convert to list type
dk_prep = dk_prep.explode(['Substrate Name','Substrate SMILES']).reset_index(drop=True)

dk_prep = dk_prep.explode(['Protein Sequence','Molecular Weight','Avg Kcat (by ec)','Avg Kcat (by ec and species)'])

dk_inp = dk_prep.drop(columns=['Gene','Reaction', 'Name', 'Avg Kcat (by ec)', 'Avg Kcat (by ec and species)','Molecular Weight'])
dk_inp.to_csv(f'../../DLKcat/DeeplearningApproach/Code/example/dk_input_{model.id}.tsv',sep="\t",na_rep='None',index= False)

logger.info("Ready for DLKcat")