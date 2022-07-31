import pandas as pd
import argparse
import re
from collections import Counter

def change_str_to_list (res_str):
    a=[]
    for i in range(len(res_str)):
        if i == 0:
            a.append(int(res_str[i][1:]))
        elif i==len(res_str)-1:
            a.append(int(res_str[i][:-1]))
        else:
            a.append(int(res_str[i]))
    return a

Metrics = ('pdb\tmain_chain\tres_num\tpfam\tpfam_name\tpdb_start\tpdb_end\tcov_num\tclan_acc\tclan_id\tuniport')
with open("./step_file/pdbbind_v2020_main_chain_PDBpfam_info.csv", 'w') as f:
    f.write(Metrics + '\n')

with open("./step_file/pdbbind_v2020_main_chain_pfam_not_find_info.csv", 'w') as f:
    f.write(Metrics + '\n')

pdb_file = pd.read_csv("/pubhome/hzhu02/GPSF/dataset/pdbbind_v2020/index/pdbbind_v2020_refine.code", header =None)
pdb_file.columns=['PDB']
pdb_list = pdb_file['PDB'].tolist()
all_chain=pd.read_csv("./step_file/pdbbind_v2020_refine_all_pocket_chain_info.csv")
mapping_file = pd.read_csv("./input_file/pdb_pfam_mapping.csv")


for pdb in pdb_list:
    sort_chain = all_chain[all_chain['pdb']==pdb].sort_values(by=['atom_num'], ascending=False)
    main_chain = sort_chain['chain'].tolist()[0]
    res_list = change_str_to_list(sort_chain['res_uniq_list'].tolist()[0].split(","))
    select_file = mapping_file[(mapping_file['PDB']==pdb)&(mapping_file['CHAIN']==main_chain)]
    
    if select_file.shape[0]>0:
        for i in range(select_file.shape[0]):
            start = select_file[i:i+1]['AUTH_PDBRES_START'].tolist()[0]
            end=select_file[i:i+1]['AUTH_PDBRES_END'].tolist()[0]
            if start != 'None' and end != 'None':
                num_start = int(select_file[i:i+1]['AUTH_PDBRES_START'].tolist()[0])
                num_end = int(select_file[i:i+1]['AUTH_PDBRES_END'].tolist()[0])
                domain = list(range(num_start, num_end+1))
                num = len(list(set(res_list).intersection(set(domain))))
                item = [pdb, main_chain, len(res_list),select_file[i:i+1]['PFAM_ACCESSION'].tolist()[0], select_file[i:i+1]['PFAM_NAME'].tolist()[0],num_start, num_end, num,select_file[i:i+1]['clan_acc'].tolist()[0],select_file[i:i+1]['clan_id'].tolist()[0],select_file[i:i+1]['UNIPROT_ACCESSION'].tolist()[0]]
            
                with open("./step_file/pdbbind_v2020_main_chain_pfam_info.csv", 'a') as f:
                    f.write("\t".join(map(str, item))+'\n')
    else:
        item = [pdb, main_chain, len(res_list),None,None,None,None,None,None,None,None ]
        with open("./step_file/pdbbind_v2020_main_chain_pfam_not_find_info.csv", 'a') as f:
            f.write("\t".join(map(str, item))+'\n')


