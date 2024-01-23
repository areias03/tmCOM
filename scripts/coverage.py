import os
import json

path = '../data/ec_data/'

for file in os.listdir(path):
    with open(os.path.join(path, file)) as f:
        ec_data = json.load(f)
        total = 0
        for k, v in ec_data.items():
            if v['mw'] == 0 or v['kcat'] == 1:
                total += 1
    print(f'File: {file}\t Coverage: {total/len(ec_data)*100}')
