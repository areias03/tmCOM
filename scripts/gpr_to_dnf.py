#!/usr/bin/env python

# Copyright (c) 2021 Moritz E. Beber
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


"""
Convert all GPR associations in a metabolic model to DNF.

: author : Moritz E. Beber
: license : BSD-3-Clause

"""


import logging
import sys
import time
import cobra
from typing import List

from cobra.io import read_sbml_model, write_sbml_model
from sympy import parse_expr
from sympy.abc import _clash1, _clash2, _clash
from sympy.logic.boolalg import to_dnf
from sympy import Basic, Float, Integer, Rational, sqrt, Symbol
from tqdm import tqdm
from reframed.io.sbml import parse_gpr_rule


logger = logging.getLogger()

local_dict= {
    "Symbol": Symbol,
    "Integer": Integer,
    "Float": Float,
    "Rational": Rational,
    "sqrt": sqrt,
}

def convert_gpr_to_dnf(model: cobra.Model) -> None:
    """
    Convert all existing GPR associations to DNF.

    RBApy can only work with the disjunctive normal form (DNF) of a
    gene-protein-reaction (GPR) association.

    """
    for rxn in tqdm(model.reactions):
        if not rxn.gene_reaction_rule:
            continue
        rule = rxn.gene_reaction_rule.replace("or", "|").replace("and", "&")
        expr = parse_expr(rule,evaluate=False,local_dict=_clash)
        rxn.gene_reaction_rule = str(to_dnf(expr)).replace("|", "or").replace("&", "and")


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
