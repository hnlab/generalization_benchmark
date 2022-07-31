import pandas as pd
import numpy as np
import oddt
from oddt import toolkit
from oddt.scoring.descriptors import close_contacts_descriptor, oddt_vina_descriptor




def extract_rf_v1_feature(ligand, protein):
    ligand_atomic_nums = [6, 7, 8, 9, 15, 16, 17, 35, 53]
    protein_atomic_nums = [6, 7, 8, 16]
    cutoff = 12
    rfscore_descriptor_cc = close_contacts_descriptor(
    cutoff=cutoff,
    protein_types=protein_atomic_nums,
    ligand_types=ligand_atomic_nums)
    feature = rfscore_descriptor_cc.build(ligand, protein)
    return rfscore_descriptor_cc.titles, feature

def extract_rf_v2_feature(ligand, protein):
    ligand_atomic_nums = [6, 7, 8, 9, 15, 16, 17, 35, 53]
    protein_atomic_nums = [6, 7, 8, 16]
    cutoff = np.array([0, 2, 4, 6, 8, 10, 12])
    descriptors = close_contacts_descriptor(
    cutoff=cutoff,
    protein_types=protein_atomic_nums,
    ligand_types=ligand_atomic_nums)
    feature = descriptors.build(ligand, protein)
    return descriptors.titles, feature


def extract_rf_vina_feature(ligand, protein):
    vina_scores = ['vina_gauss1',
                'vina_gauss2',
                'vina_repulsion',
                'vina_hydrophobic',
                'vina_hydrogen',
                'vina_num_rotors']
    vina = oddt_vina_descriptor(vina_scores=vina_scores)
    feature = vina.build(ligand, protein)
    return vina_scores, feature



pdbbind = pd.read_csv(" ", sep=",") ## load pdbid list


for i in range(pdbbind.shape[0]):
    pdb = pdbbind.iloc[i,]['pdb']
    print(pdb)
    protein = next(oddt.toolkits.ob.readfile('pdb', '/'+pdb+'/'+pdb+'_protein.pdb'))
    ligand = next(oddt.toolkits.ob.readfile('sdf', '/'+pdb+'/'+pdb+'_ligand.sdf'))
    if i == 0:
        title_1, feature_1 = extract_rf_v1_feature(ligand, protein)
        title_2, feature_2 = extract_rf_v2_feature(ligand, protein)
        title_3, feature_3 = extract_rf_vina_feature(ligand, protein)
    else:
        _, feature_36 = extract_rf_v1_feature(ligand, protein)
        _, feature_216 = extract_rf_v2_feature(ligand, protein)
        _, feature_vina = extract_rf_vina_feature(ligand, protein)
        feature_1 = np.vstack((feature_1, feature_36))
        feature_2 = np.vstack((feature_2, feature_216))
        feature_3 = np.vstack((feature_3, feature_vina))

features = np.hstack((feature_1, feature_2, feature_3))
titles = title_1 + title_2 + title_3
print(titles)
print(features)

data = pd.DataFrame(features, columns=titles)
data['pdb'] = pdbbind['pdb'].tolist()
data['affinity'] = pdbbind['affinity'].tolist()

data.to_csv("./RFScore_general_features_test.csv", index=False)
    




