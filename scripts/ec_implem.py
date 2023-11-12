#!/usr/bin/env python3

'''
Implement enzymatics constraints into models through the reaction's annotations.

:author: Alexandre Castro
'''

from reframed.io.sbml import load_cbmodel,save_cbmodel
from cobra.io import read_sbml_model, write_sbml_model
from mewpy.cobra.util import add_enzyme_constraints,convert_gpr_to_dnf,split_isozymes
from mewpy.simulation import get_simulator
from mewpy.simulation.environment import Environment as Environment
from mewpy.util.request import retreive_gene,retreive_protein,get_smiles,brenda_query

from tqdm.auto import tqdm
import pandas as pd
import numpy as np
from urllib.request import urlopen
from functools import reduce
import json
from ast import literal_eval
import requests
import argparse

parser = argparse.ArgumentParser('ec_implem')
parser.add_argument('model', help='The SBML model filepath.', type=str)
parser.add_argument('organism', help="The model's species (for better results it is not recommended to use strain name; e.g: Escherichia coli, Bacteroides fragilis).", type=str)
args = parser.parse_args()

filepath = args.model
organism = args.organism

#Reading models
model = read_sbml_model(filepath)
model2 = load_cbmodel(filepath)

#Getting sims
sim = get_simulator(model)
sim2 = get_simulator(model2)


def annotation_scraping(model_path: str):
    '''
    Run a query to find genes, reactions, substrate names, substrate SMILES and identifiers in the annotations of the model.

    :param model: The filepath to a SBML model.
    :return: A pandas.DataFrame containing the information found in the annotations.
    :rtype: pandas.DataFrame
    '''
    ls_ge = []

    #Reading models
    model = read_sbml_model(model_path)
    model2 = load_cbmodel(model_path)

    #Getting sims
    sim = get_simulator(model)
    sim2 = get_simulator(model2)


    for ge in sim.genes:
        i = sim.genes.index(ge)
        gene = sim2.genes[i]
        rxns = sim.get_gene(ge).reactions
        rx_l = []
        for rx in rxns:
            anno = sim.get_reaction(rx)['annotations']
            seed_id = anno.get('seed.reactions')
            rxn_name = sim.get_reaction(rx).name
            ecnumber = anno.get('ec-code')
            metanetx = anno.get('metanetx.reaction')
            kegg = anno.get('kegg.reaction')
            subs = list(sim.get_substrates(rx).keys())
            for sub in subs:
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
                ls_smile.append(smile)
                ls_sub.append(sub_name)
            res = [gene,rx,rxn_name,ls_sub,ls_smile,seed_id,metanetx,kegg,ecnumber]
            ls_ge.append(res)

    df_ge = pd.DataFrame(ls_ge,columns=[['Gene','Reaction','Name','Substrate Name','Substrate SMILES','ModelSEED_id','MetaNetX','KEGG_id','ecNumber']])

    return df_ge

def modelseed_query(seed_id: list):
    '''
    Queries the ModelSEED database for other id's for BIGG and KEGG.

    :param seed_id: A list or list of lists containing ModelSEED entry ID's.
    :return: List of entry ID's for BIGG, list of entry ID's for KEGG.
    :rtype: list,list
    '''
    SOLR_URL='https://modelseed.org'
    ls_kegg = []
    ls_bigg = []
    
    for mseed_id in tqdm(seed_id, 'Querying ModelSEED'):
        i = seed_id.index(mseed_id)
        if mseed_id == None:
            ls_kegg.append('None')
            ls_bigg.append('None')
        else:
            try:
                connection = urlopen(SOLR_URL+f'/solr/reactions/select?wt=json&q=id:{mseed_id}&fl=name,id,formula,charge,aliases')
                response = json.load(connection)
                for document in response['response']['docs']:
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
                        ls_bigg.append(ms_bigg)
                        ls_kegg.append(ms_kegg)
            except:
                ls_kegg.append('None')
                ls_bigg.append('None')

    return ls_bigg,ls_kegg

def bigg_query(bigg_id: list):
    '''
    Queries the BIGG database for EC numbers.
    
    :param bigg_id: A list or list of lists containing BIGG entry ID's.
    :return: List of EC numbers.
    :rtype: list
    '''
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
 
    return bigg_ls
