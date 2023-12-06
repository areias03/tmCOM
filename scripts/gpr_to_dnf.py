import logging
import sys
import time
import cobra
from typing import List

from cobra.io import read_sbml_model, write_sbml_model
from reframed.io.sbml import parse_gpr_rule, convert_to_dnf
from mewpy.simulation import get_simulator
from mewpy.util.parsing import build_tree, Boolean
from sympy import parse_expr
from sympy.logic.boolalg import to_dnf
from sympy.printing import tree
from tqdm import tqdm


logger = logging.getLogger()

def convert_gpr_to_dnf(model: cobra.Model) -> None:
    """
    Convert all existing GPR associations to DNF.
    RBApy can only work with the disjunctive normal form (DNF) of a
    gene-protein-reaction (GPR) association.
    """
    sim = get_simulator(model)
    gp = dict()
    for rxn_id in tqdm(sim.reactions):
        rxn = sim.get_reaction(rxn_id)
        if not rxn.gpr:
            continue
        gp[rxn_id] = rxn.gpr
        tree = build_tree(rxn.gpr, Boolean)
        gpr = tree.to_infix()
        print('old:',rxn.gpr)
        # TODO: update the gpr
        objective = sim.objective.get(rxn_id, 0)
        sim.add_reaction(rxn_id,
                          name=rxn.name,
                          stoichiometry=rxn.stoichiometry,
                          gpr=gpr,
                          annotations=rxn.annotations,
                          objective=objective,
                          replace=True)
        print('new:',gpr)

def convert_gpr_to_dnf2(model: cobra.Model) -> None:
    """
    Convert all existing GPR associations to DNF.
    RBApy can only work with the disjunctive normal form (DNF) of a
    gene-protein-reaction (GPR) association.
    """
    for rxn in tqdm(model.reactions):
        if not rxn.gene_reaction_rule:
            continue
        print('old:',rxn.gene_reaction_rule)
        rule = rxn.gene_reaction_rule.replace("or", "|").replace("and", "&")
        expr = build_tree(rule, Boolean)
        gpr = expr.to_infix()
        rxn.gene_reaction_rule = str(gpr).replace("|", "or").replace("&", "and")
        print('new:',rxn.gene_reaction_rule)
         
        
def main(argv: List[str]) -> None:
    """Manage the model conversion."""
    logger.info("Loading model from SBML...")
    start = time.perf_counter()
    model = read_sbml_model(sys.argv[1])
    logger.info("%.2f s", time.perf_counter() - start)
    logger.info("Converting all GPR assocations to DNF.")
    convert_gpr_to_dnf(model)
    logger.info("Writing model to SBML...")
    start = time.perf_counter()
    write_sbml_model(model, sys.argv[2])
    logger.info("%.2f s", time.perf_counter() - start)


if __name__ == "__main__":
    logging.basicConfig(level="INFO", format="[%(levelname)s] %(message)s")
    if len(sys.argv) != 3:
        logger.critical("Usage: %s <input model> <output model>", sys.argv[0])
        sys.exit(2)
    main(sys.argv[1:])
