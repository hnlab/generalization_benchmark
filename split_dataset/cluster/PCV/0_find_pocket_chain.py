import pandas as pd
import argparse
import re
from collections import Counter

def get_lig_chain(pdbid, pocket):
    z=[str(i) for i in range(0, 10)]
    pocket_file = pd.read_csv(pocket, sep="\t", header=None)
    pocket_file.columns=['all']
    pocket_file['atom']=pocket_file['all'].apply(lambda x:re.split(r"[ ]+", x)[0])
    pocket_file = pocket_file[pocket_file['atom']=="ATOM"]
    pocket_file['chain']=pocket_file['all'].apply(lambda x: re.split(r"[ ]+", x)[4]) 
    pocket_file['res']=pocket_file['all'].apply(lambda x:re.split(r"[ ]+", x)[5])
    pocket_file['new_chain'] = pocket_file['chain'].apply(lambda x: x[0])
    pocket_file['new_res']=pocket_file.apply(lambda x: x.chain[1:] if len(x.chain)>1 else x.res, axis=1)
    pocket_file['real_aa']=pocket_file['new_res'].apply(lambda x: 0 if x[-1] not in z else 1)
    chain_type = list(set(pocket_file[pocket_file['atom']=="ATOM"]['new_chain'].tolist()))

    pdb_list=[]
    chain_list=[]
    atom_list=[]
    res_num_list = []
    res_uniq_list=[]


    max_atom_num = 0
    for chain in chain_type:
        res_list = pocket_file[(pocket_file['new_chain']==chain)&(pocket_file['real_aa']==1)]['new_res'].astype(int).tolist()
        atom_num = len(res_list)
        res_list_uniq = list(set(res_list))
        res_num = len(res_list_uniq)
        
        pdb_list.append(pdbid)
        chain_list.append(chain)
        atom_list.append(atom_num)
        res_num_list.append(res_num)
        res_uniq_list.append(res_list_uniq)
    
    all_dict = {'pdb':pdb_list,'chain':chain_list,'atom_num':atom_list,'res_num':res_num_list,'res_uniq_list':res_uniq_list}
    # main_dict = {'pdb':main_pdb,'chain':main_chain,'atom_num':main_atom_num,'res_num':main_res_num,'res_uniq_list':main_uniq_res}
    all_df = pd.DataFrame(all_dict, columns=['pdb','chain','atom_num','res_num','res_uniq_list'])
    # main_df = pd.DataFrame(main_dict, columns=['pdb','chain','atom_num','res_num','res_uniq_list'])
    return all_df

def find_lig_domain(residues, select_file):

    max_num = 0
    max_index = 0
    for i in range(select_file.shape[0]):
        num_start = select_file[i:i+1]['PDB_START'].tolist()[0]
        num_end = select_file[i:i+1]['PDB_END'].tolist()[0]
        domain = list(range(num_start, num_end+1))
        num = len(list(set(residues).intersection(set(domain))))
        if num > max_num:
            max_num = num
            max_index = i
    return max_index, max_num



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="From_pdbid_to_Pfam_to_clan")
    parser.add_argument("--pdbid_list", type=str, help="path to pdbid list",
                            default="./input_file/pdbbind_v2020_general.txt") 
    parser.add_argument("--output_path",type=str, help="output path", default="./step_file")
    parser.add_argument("--clans", type=str, help="path to clans.tsv", 
                            default="./input_file/clans.tsv")
    parser.add_argument("--pdb_pfam_mapping", type=str, help="path to pdb_pfam_mapping.csv", 
                            default="./input_file/pdb_pfam_mapping.csv")
    parser.add_argument("--pocket_file", type=str, help="path to pdbid_pocket.pdb file",
                            default="./input_file")
    args = parser.parse_args().__dict__


    pdbid_file = pd.read_csv(args['pdbid_list'], header=None)
    pdbid_file.columns=['pdb']
    pdbid_list = pdbid_file['pdb'].tolist()

    pdb_chain_mapping = pd.read_csv(args["pdb_pfam_mapping"], sep=",")
    clans = pd.read_csv(args["clans"], sep="\t")
    clans = clans.rename(columns={'pfamA_acc':'PFAM_ACCESSION'})

    all_pd = pd.DataFrame(columns=['pdb','chain','atom_num','res_num','res_uniq_list'])

    for pdbid in pdbid_list:
        print(pdbid)
        pocket_file = args['pocket_file']+"/"+pdbid+"/"+pdbid+"_pocket.pdb"
        all_chain_info= get_lig_chain(pdbid, pocket_file)
        all_pd = all_pd.append(all_chain_info)


    all_pd.to_csv(args["output_path"]+"/pdbbind_v2020_general_step_1_all_pocket_chain_info.csv", index=False)


    






