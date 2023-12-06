import cobra
from typing import List

from mewpy.util.parsing import build_tree, Boolean
from tqdm import tqdm

def convert_gpr_to_dnf(model: cobra.Model) -> None:
    """
    Convert all existing GPR associations to DNF.
    RBApy can only work with the disjunctive normal form (DNF) of a
    gene-protein-reaction (GPR) association.
    """
    for rxn in tqdm(model.reactions):
        if not rxn.gene_reaction_rule:
            continue
        elif rxn.id == 'pbiosynthesis':
            continue
        print('reaction:',rxn.id)
        print('old:',rxn.gene_reaction_rule)
        rule = rxn.gene_reaction_rule.replace("or", "|").replace("and", "&")
        expr = build_tree(rule, Boolean)
        gpr = expr.to_infix()
        rxn.gene_reaction_rule = str(gpr).replace("|", "or").replace("&", "and")
        print('new:',rxn.gene_reaction_rule)        