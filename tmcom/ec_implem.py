'''
Implement enzymatics constraints into models through the reactions' annotations.

:author: Alexandre Castro
'''

from reframed.io.sbml import load_cbmodel, save_cbmodel
from cobra.io import read_sbml_model,write_sbml_model
from mewpy.cobra.util import convert_gpr_to_dnf
from mewpy.simulation import get_simulator
from mewpy.simulation.environment import Environment as Environment
from mewpy.util.request import get_smiles

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
    Run a query to find genes, reactions, substrate names, substrate SMILES and identifiers in the annotations of a SBML model.

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

def modelseed_query(seed_id: str):
    '''
    Queries the ModelSEED database for other id's for BIGG and KEGG.

    :param seed_id: A ModelSEED entry ID.
    :return: Entry ID for BIGG, list of entry ID's for KEGG.
    :rtype: str,str
    '''
    SOLR_URL='https://modelseed.org'
    
    if seed_id == None:
        bigg_id = 'None'
        kegg_id = 'None'
    else:
        try:
            connection = urlopen(SOLR_URL+f'/solr/reactions/select?wt=json&q=id:{seed_id}&fl=name,id,formula,charge,aliases')
            response = json.load(connection)
            for document in response['response']['docs']:
                bigg_id = list(filter(lambda a: 'BiGG:' in a, document.get('aliases')))
                kegg_id = list(filter(lambda a: 'KEGG:' in a, document.get('aliases')))
                if len(bigg_id)== 0 and len(kegg_id)== 0:
                    bigg_id = 'None'
                    kegg_id = 'None'
                elif len(bigg_id)== 0 and len(kegg_id)!= 0:
                    bigg_id = 'None'
                    kegg_id = list(kegg_id)[0]
                    kegg_id = kegg_id.replace('KEGG: ','')
                elif len(bigg_id)!= 0 and len(kegg_id)== 0:
                    kegg_id = 'None'
                    bigg_id = list(bigg_id)[0]
                    bigg_id = bigg_id.replace('BiGG: ','')
                else:
                    kegg_id = list(kegg_id)[0]
                    kegg_id = kegg_id.replace('KEGG: ','')
                    bigg_id = list(bigg_id)[0]
                    bigg_id = bigg_id.replace('BiGG: ','')
        except:
            bigg_id = 'None'
            kegg_id = 'None'

    return bigg_id,kegg_id

def bigg_query(bigg: str):
    '''
    Queries the BIGG database for EC numbers.
    
    :param bigg: A BIGG entry ID.
    :return: List of EC numbers.
    :rtype: list
    '''
    
    bigg_ls = []
    url =f'http://bigg.ucsd.edu/api/v2/universal/reactions/{bigg}'

    with requests.request("GET", url) as resp:
        try:
            resp.raise_for_status()  # raises exception when not a 2xx response
            if resp.status_code != 204:
                data = dict(resp.json())
                ec_l = data['database_links']
                if ec_l == None:
                    bigg_ls.append(None)
                else:
                    ec = [i['id'] for i in ec_l['EC Number']]
                    if ec == None:
                        bigg_ls.append(None)
                    bigg_ls.append(ec)
            else: 
                bigg_ls.append(None)
        except:
            bigg_ls.append(None)

    return bigg_ls
