from mewpy.model.com import CommunityModel
from cobra.io import read_sbml_model
from mewpy.cobra.com import SteadyCom
from concurrent.futures import ProcessPoolExecutor


ec_bt = read_sbml_model('../models/ec/ec_Bacteroides_thetaiotaomicron_VPI_5482.xml')
ec_ec = read_sbml_model('../models/ec/ec_Escherichia_coli_ED1a.xml')

sample = CommunityModel([ec_bt, ec_ec], abundances=[1, 0.5])


def main(val: int):
    constraints = {'R_EX_fe2[e]': (-0.00007, 0), 'R_EX_pheme[e]': (-0.0000001, 0)
                   , 'R_prot_pool_exchange_M_Bacteroides_thetaiotaomicron_VPI_5482': (0, val)
                   , 'R_prot_pool_exchange_M_Escherichia_coli_ED1a': (0, 0.095)}
    try:
        SteadyCom(sample, constraints=constraints)
        print(f'Grew at {val}')
    except ZeroDivisionError:
        print(f'Did not grow at {val}')


if __name__ == "__main__":
    with ProcessPoolExecutor() as p:
        p.map(main, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
