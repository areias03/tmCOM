from mewpy.model.com import CommunityModel
from cobra.io import read_sbml_model
from tmcom.steadycom import SteadyCom
from concurrent.futures import ProcessPoolExecutor
import itertools
from typing import List, Tuple

ec_bt = read_sbml_model('models/ec/ec_Bacteroides_thetaiotaomicron_VPI_5482.xml')
ec_ec = read_sbml_model('models/ec/ec_Escherichia_coli_ED1a.xml')

sample = CommunityModel([ec_bt, ec_ec])


def main(args: Tuple[float, float]):
    val_1, val_2 = args
    constraints = {'R_EX_fe2[e]': (-0.00007, 0), 'R_EX_pheme[e]': (-0.0000001, 0)
                   , 'R_prot_pool_exchange_M_Bacteroides_thetaiotaomicron_VPI_5482': (0, val_2)
                   , 'R_prot_pool_exchange_M_Escherichia_coli_ED1a': (0, val_1)}
    try:
        sol = SteadyCom(sample, constraints=constraints)
        res  = f'Grew at Ec:{val_1}, Bt:{val_2}\n {sol}\n'
    except ZeroDivisionError:
        res = f'Did not grow at Ec:{val_1}, Bt:{val_2}\n'
    return res


if __name__ == "__main__":

    val_1 = [0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3]
    val_2 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    args: List[Tuple[float,float]] = list(itertools.product(val_1, val_2))
    print('Possible combinations:', f"Value 1: {val_1} | Value 2: {val_2} ")
    with ProcessPoolExecutor() as p:
        with open('data/results/prot_pool.txt','w') as f:
            for result in p.map(main, args):
                print(result, file = f)
