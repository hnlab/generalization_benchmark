from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA
import mdtraj as md
import numpy as np
import pandas as pd
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_absolute_error
from scipy.stats import spearmanr, pearsonr
import seaborn as sns
from scipy.stats import gaussian_kde
import oddt
from sklearn.ensemble import RandomForestRegressor
import pickle
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import PredefinedSplit
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from oddt.scoring import scorer, ensemble_model
from oddt.utils import method_caller
from oddt.scoring.models.regressors import neuralnetwork
import xgboost as xgb
from sklearn.linear_model import SGDRegressor
from sklearn.inspection import permutation_importance
import shap
from collections import defaultdict
# import pandas as pd
# import numpy as np
import oddt
from oddt import toolkit
from oddt.scoring.descriptors import close_contacts_descriptor, oddt_vina_descriptor

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from sklearn.datasets import load_breast_cancer
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
import sys, argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='calculate buried SASA of complex structures')
    parser.add_argument('--input',type=int)
    parser.add_argument('--output',type=str)
    args = parser.parse_args().__dict__


    all_pdb = pd.read_csv("../PCV/input_file/pdbbind_v2020_general.txt", header=None)
    all_pdb.columns=['pdb']

    Metrics = ('pdb\tligand_sasa\tprotein_sasa\tcomplex_sasa\tdelta_sasa')

    pdb_list = all_pdb['pdb'].to_list()
    with open(args["output"], 'w') as f:
        f.write(Metrics + '\n')

    for pdb in pdb_list[args["input"]:args["input"]+1]:

        try:
            ligand = md.load("./"+pdb+"/"+pdb+"_ligand_ob_new.pdb")
            protein = md.load("./"+pdb+"/"+pdb+"_protein_new.pdb")
            complex = md.load("./"+pdb+"/"+pdb+"_complex_new.pdb")
            lig_sa = np.sum(md.shrake_rupley(ligand))
            pro_sa = np.sum(md.shrake_rupley(protein))
            com_sa = np.sum(md.shrake_rupley(complex))
            del_sa = lig_sa+pro_sa - com_sa
            with open(args["output"], 'a') as f:
                f.write("\t".join(map(str, [pdb,lig_sa,pro_sa,com_sa,del_sa]))+'\n')
        except:
            print(pdb)

