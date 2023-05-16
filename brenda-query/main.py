import pandas as pd
from mewpy.util.request import brenda_query


df = pd.read_csv('./rxn_data.csv')

kcat_ls = []

for ec in df['ecNumber'].values.tolist():
    sub_kcat_ls = []
    for i in range(len(ec)):
        ec_n = ec[i]
        kcat = kcat = brenda_query(user = 'pg45962@uminho.pt',password='Mentafrio+15',ecNumber=ec_n,organism = 'Escherichia coli')
        sub_kcat_ls.append(kcat)
    kcat_ls.append(sub_kcat_ls)                
