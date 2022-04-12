import pandas as pd

file = pd.read_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/pfam/general/pdbbind_v2020_general_step_1_all_pocket_chain_info.csv")
general_code = pd.read_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/pfam/general/pdbbind_v2020_general.txt", header=None)
general_code.columns=['pdb']

main_chain = pd.DataFrame()
for code in general_code['pdb'].tolist():
    subfile = file[file['pdb']==code]
    subfile.sort_values(by=['atom_num'], ascending=False)
    main_chain = pd.concat([main_chain, subfile.head(1)])

main_chain.columns=file.columns.tolist()

main_chain.to_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/pfam/general/pdbbind_v2020_general_main_chain.csv", index=False)

