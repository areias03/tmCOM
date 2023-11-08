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


def annotation_scraping(model_path: Union["Model", "CBModel"]):
    '''
    Run a query to find gens, reactions and identifiers in the annotations of the model.

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
            res = [gene,rx,rxn_name,seed_id,metanetx,kegg,ecnumber]
            ls_ge.append(res)

    df_ge = pd.DataFrame(ls_ge,columns=[['Gene','Reaction','Name','ModelSEED_id','MetaNetX','KEGG_id','ecNumber']])

    return df_ge

def modelseed_query(seed_id):
    '''
    Queries the ModelSEED database for other id's for BIGG and KEGG.

    :param df_ge: A pandas.DataFrame containing the information found in the annotations.
    :return: The input pandas.DataFrame with adicional entry ID's to other databases (BIGG and KEGG), the list of ID's for BIGG.
    :rtype: pandas.DataFrame
    '''
    #Converting ID's to list
    seed_id = df_ge['ModelSEED_id'].values.tolist()
    seed_id = [reduce(lambda x: x, inner_list) for inner_list in seed_id]

    kegg_id = df_ge['KEGG_id'].values.tolist()
    kegg_id = [reduce(lambda x: x, inner_list) for inner_list in kegg_id]

    #Querying DB
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

    return df_ge,ls_bigg
