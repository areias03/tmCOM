import argparse
import os
import sys
from typing import List

import numpy as np
from cobra.io import read_sbml_model
from cobra.io.sbml import Model
from mewpy.simulation import Environment, Simulator, get_simulator


def media_growth_test(model: Model):
    print(model)


if __name__ == "__main__":
    pass
