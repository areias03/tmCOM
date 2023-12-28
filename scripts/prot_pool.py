import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec

path = '../data/model_tuning/'

files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
files = sorted(files)
ind_list = []
merged = []

l_va = np.linspace(0.1, 1.0, num=10)

for i in range(len(l_va)):
        ind_list.append(f'{round(l_va[i], 1)}')
        
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
df = df.loc[:,~df.columns.duplicated()]
print(df)
'''

# Generate clustermap

x_labels = ['Bt VPI-5482','Bu ATCC-8492','Ec ED1a','Fn ATCC-25586','Ri L1-82','Sp ATCC-15912','Ss DSM-20560']
g = sns.clustermap(df, row_cluster = False, vmin = 0, xticklabels = x_labels)
g.tick_params(axis='both', which='both', length=0)
g.ax_heatmap.set_ylabel(g.ax_heatmap.get_ylabel(), labelpad = 20)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)  # For y axis
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90) # For x axis
plt.savefig('../data/images/prot_pool_all.png', bbox_inches='tight')
'''

df_inv = df.T
print(df_inv)

# Generate violin plot
v = sns.violinplot(data=df_inv, inner="point", palette="Blues")
v.tick_params(axis='both', which='both', length=0)
v.set_ylabel("Growth", labelpad = 15)
v.set_xlabel(v.get_xlabel(), labelpad = 15)
plt.setp(v.get_yticklabels(), rotation=0)  # For y axis
plt.setp(v.get_xticklabels(), rotation=90) # For x axis
plt.savefig('../data/images/prot_pool_violin.png', bbox_inches='tight')

'''
df_inv = df.T
print(df_inv)

x_labels = ['Bt VPI-5482','Bu ATCC-8492','Ec ED1a','Fn ATCC-25586','Ri L1-82','Sp ATCC-15912','Ss DSM-20560']

#First create the clustermap figure
g = sns.clustermap(df, row_cluster = False, vmin = 0, xticklabels = x_labels, figsize=(15,10))
g.tick_params(axis='both', which='both', length=0)
g.ax_heatmap.set_ylabel(g.ax_heatmap.get_ylabel(), labelpad = 20)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)  # For y axis
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90) # For x axis
# set the gridspec to only cover half of the figure
g.gs.update(left=0.05, right=0.45)

#create new gridspec for the right part
gs2 = matplotlib.gridspec.GridSpec(1,1, left=0.6)
# create axes within this new gridspec
ax2 = g.fig.add_subplot(gs2[0])
# plot boxplot in the new axes
v = sns.violinplot(data=df_inv, inner="point", ax = ax2)
v.tick_params(axis='both', which='both', length=0)
v.set_ylabel("Growth", labelpad = 15)
v.set_xlabel(v.get_xlabel(), labelpad = 15)
plt.setp(v.get_yticklabels(), rotation=0)  # For y axis
plt.setp(v.get_xticklabels(), rotation=90) # For x axis
plt.savefig('../data/images/prot_pool.png', bbox_inches='tight')
'''