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

def convert_to_irreversible(model: Union[Simulator, "Model", "CBModel"], inline: bool = False):
    """Split reversible reactions into two irreversible reactions
    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    :param model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator
    :return: a irreversible model simulator, a reverse mapping.
    :rtype:(Simulator,dict)
    """
    
    sim = get_simulator(deepcopy(model))

    objective = sim.objective.copy()
    irrev_map=dict()
    for r_id in tqdm(sim.reactions, "Converting to irreversible"):
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            rxn = sim.get_reaction(r_id)
            rev_rxn_id = r_id+"_REV"

            rev_rxn = dict()
            rev_rxn['name'] = rxn.name + " reverse"
            rev_rxn['lb'] = 0
            rev_rxn['ub'] = -rxn.lb
            rev_rxn['gpr'] = rxn.gpr
            sth = {k: v * -1 for k, v in rxn.stoichiometry.items()}
            rev_rxn['stoichiometry'] = sth
            rev_rxn['reversible'] = False
            rev_rxn['annotations'] = copy(rxn.annotations)

            sim.add_reaction(rev_rxn_id, **rev_rxn)
            sim.set_reaction_bounds(r_id, 0, rxn.ub, False)
            
            irrev_map[r_id] = rev_rxn_id
            
            if r_id in objective:
                objective[rev_rxn_id] = -objective[r_id]

    sim.objective = objective
    return sim, irrev_map


def split_isozymes(model: Union[Simulator, "Model", "CBModel"], inline: bool = False):
    """Splits reactions with isozymes into separated reactions

    :param model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator
    :param (boolean) inline: apply the modifications to the same of generate a new model. Default generates a new model.
    :return: a simulator and a mapping from original to splitted reactions
    :rtype: (Simulator, dict)
    """

    sim = get_simulator(deepcopy(model))

    objective = sim.objective
    mapping = dict()
    newobjective = {}

    for r_id in tqdm(sim.reactions, "Splitting isozymes"):
        rxn = sim.get_reaction(r_id)
        gpr = rxn.gpr

        if gpr is not None and len(gpr.strip()) > 0:
            t = build_tree(rxn.gpr, Boolean)
            ren = {x: r_gene(x) for x in t.get_operands()}
            gpr = t.replace(ren).to_infix()
            proteins = isozymes(gpr)
            mapping[r_id] = []
            for i, protein in enumerate(proteins):
                r_id_new = '{}_No{}'.format(r_id, i+1)
                mapping[r_id].append(r_id_new)

                rxn_new = dict()
                rxn_new['name'] = '{} No{}'.format(rxn.name, i+1)
                rxn_new['lb'] = rxn.lb
                rxn_new['ub'] = rxn.ub
                rxn_new['gpr'] = protein
                rxn_new['stoichiometry'] = rxn.stoichiometry.copy()
                rxn_new['annotations'] = copy(rxn.annotations)

                sim.add_reaction(r_id_new, **rxn_new)
            sim.remove_reaction(r_id)

    # set the objective
    for k, v in objective.items():
        if k in mapping.keys():
            for r in mapping[k]:
                newobjective[r] = v
        else:
            newobjective[k] = v

    sim.objective = newobjective

    return sim, mapping


def __enzime_constraints(model: Union[Simulator, "Model", "CBModel"],
                         prot_mw=None,
                         enz_kcats=None,
                         c_compartment: str = 'c',
                         inline: bool = False):
    """Auxiliary method to add enzyme constraints to a model

    :param model: A model or simulator
    :type model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator
    :param data: Protein MW and Kcats 
    :type data: None
    :param c_compartment: The compartment where gene/proteins pseudo species are to be added.
        Defaults to 'c'
    :type c_compartment: str, optional
    :param (boolean) inline: apply the modifications to the same of generate a new model.
        Default generates a new model.
    :type inline: bool, optional
    :return: a new enzyme constrained model
    :rtype: Simulator
    """

    if inline:
        sim = get_simulator(model)
    else:
        sim = deepcopy(get_simulator(model))
    
    objective = sim.objective

    if prot_mw is None:
        prot_mw = dict()
        for gene in sim.genes:
            prot_mw[gene] = {'protein': gene[len(sim._g_prefix):], 'mw': 1}

    if enz_kcats is None:
        enz_kcats = dict()
        for gene in sim.genes:
            enz_kcats[gene] = dict()
            rxns = sim.get_gene(gene).reactions
            for rxn in rxns:
                enz_kcats[gene][rxn] = {'protein': gene[len(sim._g_prefix):], 'kcat': 1}


    # Add protein pool and species
    common_protein_pool_id = sim._m_prefix+'prot_pool_c'
    pool_reaction = sim._r_prefix+'prot_pool_exchange'

    sim.add_metabolite(common_protein_pool_id,
                       name='prot_pool [cytoplasm]',
                       compartment=c_compartment)

    sim.add_reaction(pool_reaction,
                     name='protein pool exchange',
                     stoichiometry={common_protein_pool_id: 1},
                     lb=0,
                     ub=inf,
                     reversible=False,
                     reaction_type='EX'
                     )

    # Add gene/protein species and draw protein pseudo-reactions
    # MW in kDa, [kDa = g/mmol]
    gene_meta = dict()
    for gene in tqdm(sim.genes, "Adding gene species"):
        mw = prot_mw[gene]
        m_prot_id = f"prot_{mw['protein']}_{c_compartment}"
        m_name = f"prot_{mw['protein']} {c_compartment}"
        sim.add_metabolite(m_prot_id,
                           name=m_name,
                           compartment=c_compartment)

        gene_meta[gene] = m_prot_id

        r_prot_id = f"draw_prot_{mw['protein']}"
        sim.add_reaction(r_prot_id,
                         name=r_prot_id,
                         stoichiometry={common_protein_pool_id: -1*mw['mw'],
                                        m_prot_id: 1},
                         lb=0,
                         ub=inf,
                         reversible=False,
                         )

    # Add enzymes to reactions stoichiometry.
    # 1/Kcats in per hour. Considering kcats in per second.
    for rxn_id in tqdm(sim.reactions, "Adding proteins usage to reactions"):
        rxn = sim.get_reaction(rxn_id)
        if rxn.gpr:
            s = rxn.stoichiometry
            genes = build_tree(rxn.gpr, Boolean).get_operands()
            for g in genes:
                # TODO: mapping of (gene, reaction ec) to kcat
                try:
                    s[gene_meta[g]] = -1/(enz_kcats[g][rxn_id]['kcat']*3600)
                except Exception:
                    s[gene_meta[g]] = -1/(ModelConstants.DEFAULT_KCAT*3600)
            sim.update_stoichiometry(rxn_id, s)
    sim.objective = objective
    return sim


def add_enzyme_constraints(model: Union[Simulator, "Model", "CBModel"],
                           prot_mw=None,
                           enz_kcats=None,
                           c_compartment: str='c',
                           inline: bool=False):
    """Adds enzyme constraints to a model.

    :param model: A model or simulator
    :type model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator
    :param data: Protein MW and Kcats 
    :type data: None
    :param c_compartment: The compartment where gene/proteins pseudo species are to be added.
        Defaults to 'c'
    :type c_compartment: str, optional
    :param (boolean) inline: apply the modifications to the same of generate a new model.
        Default generates a new model.
    :type inline: bool, optional
    :return: a new enzyme constrained model
    :rtype: Simulator
    """
    sim, _ = convert_to_irreversible(model, inline)
    sim, _ = split_isozymes(sim, True)
    sim = __enzime_constraints(sim,
                               prot_mw=prot_mw,
                               enz_kcats=enz_kcats,
                               c_compartment=c_compartment,
                               inline=True)
    return sim
