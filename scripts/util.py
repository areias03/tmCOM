from reframed.io.sbml import load_cbmodel
from cobra.io import read_sbml_model
from mewpy.simulation import get_simulator
from mewpy.util.request import get_smiles,brenda_query
from tqdm.auto import tqdm
import pandas as pd
import numpy as np
from urllib.request import urlopen
from functools import reduce
import json
from ast import literal_eval
import requests
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from cobra import Model
    from reframed.core.cbmodel import CBModel


def annotation_scraping(model_path: str):
    '''
    Run a query to find genes, reactions, substrate names, substrate SMILES and identifiers in the annotations of the model.

    :param model: The filepath to a SBML model.
    :type model_path: str
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


    for ge in tqdm(sim.genes):
        i = sim.genes.index(ge)
        gene = sim2.genes[i]
        rxns = sim.get_gene(ge).reactions
        for rx in rxns:
            anno = sim.get_reaction(rx)['annotations']
            seed_id = anno.get('seed.reactions')
            rxn_name = sim.get_reaction(rx).name
            sub = list(sim.get_substrates(rx).keys())
            ecnumber = anno.get('ec-code')
            metanetx = anno.get('metanetx.reaction')
            kegg = anno.get('kegg.reaction')
            res = [gene,rx,rxn_name,seed_id,metanetx,kegg,ecnumber]
            ls_ge.append(res)

    df_ge = pd.DataFrame(ls_ge,columns=[['Gene','Reaction','Name','Substrates','ModelSEED_id','MetaNetX','KEGG_id','ecNumber']])

    return df_ge

def modelseed_query(seed_id: list):
    '''
    Queries the ModelSEED database for other id's for BIGG and KEGG.

    :param seed_id: A list or list of lists containing ModelSEED entry ID's.
    :type seed_id: list
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
    :type biig_id: list
    :return: List of EC numbers.
    :rtype: list
    '''
    bigg_ls = []

    for bigg in tqdm(ls_bigg, 'Querying BIGG'):
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

def find_subs_smiles(model: Union["Model", "CBModel"]):
    '''
    Finds substrate ID's, names and canonical SMILES, and reviews substrate names to find more SMILES.
    
    :param model: A model.
    :type model: A COBRApy or REFRAMED Model.
    :return: List of substrates, list of names, list of SMILES.
    :rtype: list,list,list
    '''
  
    sim = get_simulator(model)
    ls_sub_n = []
    ls_smile =[]
    ls_sub = []
    
    for ge in tqdm(sim.genes):
        rxns = sim.get_gene(ge).reactions
        for rx in rxns:
            subs = list(sim.get_substrates(rx).keys())
            ls_sub.append(subs)
            sub_ls_sub = []
            sub_ls_smile = []
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
                sub_ls_smile.append(smile)
                sub_ls_sub.append(sub_name)
            ls_sub_n.append(sub_ls_sub)
            ls_smile.append(sub_ls_smile)
                
    return ls_sub,ls_sub_n,ls_smile