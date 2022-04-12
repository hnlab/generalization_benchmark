import pandas as pd
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
import re

def get_chain(pdb_id, chain):
    with open("/pubhome/hzhu02/GPSF/dataset/pdbbind_v2020/v2020-other-PL/"+pdb_id+"/"+pdb_id+"_protein.pdb", "r") as f:
        chain_seq = {}
        for line in f.readlines():
            line = line.strip().split()
            if line[0]=="ATOM" and line[4][0]==chain:
                if len(line[4]) == 1:
                    num = re.sub("[A-Z]","",line[5])
                else:
                    num = line[4][1:]
                    
                if num not in chain_seq.keys():
                    aa =line[3]
                    if is_aa(aa, standard=True):
                        one_aa = three_to_one(aa)
                    elif aa in {'HIE','HID'}:
                        one_aa = 'H'
       

                    elif aa in {'CYX','CYM'}:
                            one_aa = 'C'
                    else:
                        one_aa = 'X'

                    chain_seq[num]=one_aa
        
    return chain_seq


not_found_pdb = pd.read_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/pfam/general/not_find_pfam.code", header=None)
not_found_pdb.columns=['pdb']

pdb_list = not_found_pdb['pdb'].tolist()
all_chain = pd.read_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/pfam/general/pdbbind_v2020_general_step_1_all_pocket_chain_info.csv")

num=0
for pdb in pdb_list:
    subfile = all_chain[all_chain['pdb']==pdb]
    subfile = subfile.sort_values(by=['atom_num'], ascending=False)
    chain = subfile['chain'].tolist()[0]
    # print(pdb,chain)

    chain_seq = get_chain(pdb, chain)
    # start = int(min(chain_seq.keys()))
    # print(start)
    # end = int(max(chain_seq.keys()))
    # print(end)
    # counts = [i for i in range(pdb_start, pdb_end+1)]
    seq=""

    for i in chain_seq.keys():
        try:
            # print(i)
            seq += chain_seq[i]
        except:
            # print("worning pdn seq not continue", pdb)
            num += 1

    print(">"+pdb)
    print(seq)
 



            
