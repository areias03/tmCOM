import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

path = '../data/model_tuning/'

files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
ind_list = []
merged = []

l_va = np.linspace(0.1, 1.0, num=10)

for i in range(len(l_va)):
        ind_list.append(f'0 - {round(l_va[i], 1)}')
        
ind = {'Protein pool exchange': ind_list}

df_ind = pd.DataFrame(ind)

merged.append(df_ind)

for p in files:
    p = path + p
    print(p)
    r = pd.read_csv(p)
    r = r.drop(r.columns[0],axis=1)
    merged.append(r)
df = pd.concat(merged, axis=1)
df.to_csv('../data/model_tuning/results/prot_pool_all.csv')
df = df.set_index('Protein pool exchange')
print(df)
#df = df.drop(df.columns[0],axis=1)
#df_norm = df.apply(lambda x: (x-x.mean())/x.std(), axis = 1)
#print(df_norm)
sns.heatmap(df)
plt.savefig('../data/model_tuning/results/prot_pool_all.png', bbox_inches='tight')
    
'''
for result in p.map(prot_pool_analysis, m_list):
            merged.append(result)
        df = pd.concat(merged, axis=1)
        df.to_csv('../data/model_tuning/prot_pool_all.csv')
        print(df)
        df_norm = (df - df.mean()) / df.std()
        sns.heatmap(df_norm)
        plt.savefig('../data/model_tuning/prot_pool_all.png')
'''