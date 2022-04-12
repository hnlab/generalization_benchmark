import pandas as pd

all_pfam = pd.read_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/pfam/general/pdbbind_v2020_main_chain_PDBPfam_info.csv", sep="\t")

pdb_file = pd.read_csv("/pubhome/hzhu02/GPSF/dataset/pdbbind_v2020/pdbbind_v2020_general.txt", header =None)
pdb_file.columns=['PDB']
pdb_list = pdb_file['PDB'].tolist()
# print(all_pfam.columns)

results=pd.DataFrame()
# to_be_check = pd.DataFrame()
not_found=[]
for pdb in pdb_list:
    print(pdb)
    pfams = all_pfam[(all_pfam['pdb']==pdb)&(all_pfam['cov_num']!=0)]
    if pfams.shape[0] ==1:
        results=pd.concat([results, pfams])
    elif pfams.shape[0] >1 :
        pfams = pfams.sort_values(by=['cov_num'], ascending=False)
        values = pfams['cov_num'].tolist()
        res_num = pfams['res_num'].tolist()[0]
        persent = float(float(values[0])/float(res_num))
        print(persent)

        results = pd.concat([results, pfams.head(1)])
    else:
        not_found.append(pdb)

results.to_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/pfam/general/pdbbind_v2020_main_pfam.csv", index=False)
# to_be_check.to_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/pfam/general/pdbbind_v2020_main_pfam_to_be_check.csv", index=False)
not_found = pd.DataFrame(not_found)
not_found.to_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/pfam/general/pdbbind_v2020_pfam_not_found.csv", index=False)