import pandas as pd
result = []
cluster_type=[]
with open("/pubhome/hzhu02/Redocked_pose/split_dataset/3_fold_ccv/general_refine/general_refine_cluster.out", "r") as f:
    for line in f.readlines():
        line = line.strip().split(":")
        cluster = line[0]
        cluster_type.append(cluster)
        # print(cluster)
        pdbs = line[1].strip().split(" ")
        for pdb in pdbs:
            # print(pdb)
            result.append([cluster, pdb])

result = pd.DataFrame(result, columns=['cluster','pdb'])
result.to_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/3_fold_ccv/general_refine/SCV_general_refine.csv", index=False)

summary_result = []
for item in cluster_type:
    item_file = result[result['cluster']==item]
    summary_result.append([item, item_file.shape[0]])

summary_result = pd.DataFrame(summary_result, columns=['cluster','pdb_num'])
summary_result = summary_result.sort_values(by=['pdb_num'])
summary_result.to_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/3_fold_ccv/general_refine/general_refine_SCV_summary.csv", index=False)
