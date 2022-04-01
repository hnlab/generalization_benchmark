#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ========================================================================
# partiton the dataset using cluster_cross_validation
# modified from gnina   https://github.com/gnina/scripts/blob/master/clustering.py
# ========================================================================

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio import pairwise2
from multiprocessing import Pool, cpu_count
from functools import partial
from matplotlib.pyplot import get
import scipy.cluster.hierarchy
import numpy as np
import pandas as pd
import sys, argparse, bisect, re, os, fnmatch
import pickle, collections
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import FingerprintSimilarity as fs
from rdkit.Chem.Fingerprints import FingerprintMols
import rdkit

''' Some amino acids have nonstandard residue names: 
http://ambermd.org/tutorials/advanced/tutorial1_adv/
HIE histidine H
HID histidine H
HEM hemoglobin? (http://www.bmsc.washington.edu/CrystaLinks/man/pdb/guide2.2_frame.html)
CYX cystenine C
CYM cystenine C'''

def computeLigandSimilarity(pdbid_list, fname):
    '''Read target (first col) and ligand (third col) from fname. 
    Return ligand similarity matrix indexed according to pdbid_list'''
    fingerprints = dict()
    for line in open(fname):
        vals = line.split()
        targ = vals[0]
        ligfile = vals[2]
        smi = open(ligfile).readline().split()[0]
        mol = AllChem.MolFromSmiles(smi)
        if mol == None:
            mol = AllChem.MolFromSmiles(smi,sanitize=False)
        fp = FingerprintMols.FingerprintMol(mol)
        fingerprints[targ] = fp
    n = len(pdbid_list)
    sims = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1):
            fpi = fingerprints[pdbid_list[i]]
            fpj = fingerprints[pdbid_list[j]]
            sim = fs(fpi,fpj)
            sims[i,j] = sims[j,i] = sim
    return sims

def getResidueStrings(structure):
    seqs = []
    for model in structure:
        for ch in model.get_chains():
            seq = ''
            for residue in model.get_residues():
                resname = residue.get_resname()
                if is_aa(resname, standard=True):
                    seq += three_to_one(resname)
                elif resname in {'HIE', 'HID'}:
                    seq += 'H'
                elif resname in {'CYX', 'CYM'}:
                    seq += 'C'
                else:
                    seq += 'X'
            seqs.append(seq)
    return seqs

def cUTDM2(targets, pair):
    '''compute distance between target pair'''
    (a, b) = pair
    mindist = 1.0
    for seq1 in targets[a]:
        for seq2 in targets[b]:
            score = pairwise2.align.globalxx(seq1, seq2, score_only=True)
            length = max(len(seq1), len(seq2))
            distance = (length-score)/length
            if distance < mindist:
                mindist = distance
    #print (a,b,mindist)
    return (a, b, mindist)

def calProteinDistanceMatrix(targets):
    '''compute full pairwise target distance matrix in parallel'''
    n = len(targets)
    pairs = [(r, c) for r in range(n) for c in range(r+1, n)] #upper triangle
    pool = Pool()
    function = partial(cUTDM2, targets)
    distanceTuples = pool.map(function, pairs)
    distanceMatrix = np.zeros((n, n))
    for (a, b, distance) in distanceTuples:
        distanceMatrix[a][b] = distanceMatrix[b][a] = distance
    return distanceMatrix

def calLigandSimilarity(fingerprints, pdbid):
    n = len(pdbid)
    sims = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1):
            fpi = fingerprints[pdbid[i]]
            fpj = fingerprints[pdbid[j]]
            sim = fs(fpi,fpj)
            sims[i,j] = sims[j,i] = sim
    return sims


def assignGroup(dists, ligandsim, t, t2, ligandt, explore, names):
    '''group targets that are less than t away from each other and what's in explore'''
    group = set(explore)
    while explore:
        frontier = set()
        for i in explore:
            for j in range(dists.shape[1]):
                if j not in group:
                    #add to the group if protein is close by threshold t (these are distances - default 0.5)
                    #also add if the ligands are more similar (not distance) than ligandt and 
                    #the protein is closer than t2 (default 0.8 - meaning more than 20% similar)
                    if dists[i][j] < t or (ligandsim[i][j] > ligandt and dists[i][j] < t2):
                        group.add(j)
                        frontier.add(j)                
                                        
        explore = frontier
    return group

def calcClusterGroups(dists, ligandsim, pdbid_list, t, t2, ligandt):
    '''dists is a distance matrix (full) for pdbid_list'''
    assigned = set()
    groups = []
    for i in range(dists.shape[0]):
        if i not in assigned:
            group = assignGroup(dists, ligandsim, t, t2, ligandt, set([i]),pdbid_list)
            groups.append(group)
            assigned.update(group)
    return [set(pdbid_list[i] for i in g) for g in groups]

def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()

def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    else: return -1

def createFolds(cluster_groups, numfolds, randomize):
    '''split target clusters into numfolds folds with balanced num poses per fold
       If randomize, will balance less well.
    '''
    folds = [[] for _ in range(numfolds)]
    fold_numposes = [0]*numfolds
    group_numposes = [0]*len(cluster_groups)
    foldmap = {}
    for i, group in enumerate(cluster_groups):
        for target in group:            
            group_numposes[i] += 1
    for _ in cluster_groups:
        #iteratively assign group with most poses to fold with fewest poses
        maxgroup = group_numposes.index(np.max(group_numposes))
        if randomize:
            space = np.max(fold_numposes) - np.array(fold_numposes)
            tot = np.sum(space)
            if tot == 0:
                minfold = np.random.choice(numfolds)
            else: #weighted selection, prefer spots with more free space
                choice = np.random.choice(tot)
                tot = 0
                for i in range(len(space)):
                    tot += space[i]
                    if choice < tot:
                        minfold = i
                        break
        else:
            minfold = fold_numposes.index(np.min(fold_numposes))
        folds[minfold].extend(cluster_groups[maxgroup])
        fold_numposes[minfold] += group_numposes[maxgroup]
        group_numposes[maxgroup] = -1
        for t in cluster_groups[maxgroup]:
            foldmap[t] = minfold
    print('Poses per fold: {}'.format(fold_numposes))
    for f in folds:
        f.sort()
    return folds, foldmap


def checkFolds(dists, target_names, threshold, foldmap):
    '''check that targets in different folds pass dissimilarity threshold'''
    ok = True
    n_targets = dists.shape[0]
    min_dist = np.inf
    closest = None
    for t in foldmap:
        if t not in set(target_names):
            print('warning: {} not found in distance matrix'.format(t))
    for a in range(n_targets):
        for b in range(a+1, n_targets):
            a_name = target_names[a]
            b_name = target_names[b]
            if a_name in foldmap and b_name in foldmap:
                if foldmap[a_name] != foldmap[b_name]:
                    if dists[a][b] < min_dist:
                        min_dist = dists[a][b]
                        closest = (a_name, b_name)
                    if dists[a][b] < threshold:
                        print('warning: {} and {} are {:.3f}% similar but in different folds' \
                              .format(a_name, b_name, 100*(1-dists[a][b])))
                        ok = False
    if closest:
        print('{} and {} are the most similar targets in different folds ({:.3f}%)' \
              .format(closest[0], closest[1], 100*(1-min_dist)))
    return ok

def loadFolds(inname, target_names, numfolds):
    #load test/train files
    trainfiles = [open('{}train{}.types'.format(inname,x),'r') for x in range(numfolds)]
    testfiles = [open('{}test{}.types'.format(inname,x),'r') for x in range(numfolds)]
    folds = [set() for _ in range(numfolds)]
    foldmap = {}
    target_set = set(target_names)
    for i in range(numfolds):
        for line in testfiles[i]:
            for word in re.findall(r'[\w]+', line):
                if word in target_set:
                    target = word
                    break
            for j in range(numfolds):
                if j == i:
                    folds[i].add(target)
                else:
                    assert target not in folds[j]
            foldmap[target] = i
        for line in trainfiles[i]:
            for word in re.findall(r'[\w]+', line):
                if word in target_set:
                    target = word
                    break
            assert target not in folds[i]
    return folds, foldmap

def get_chain(pdb_id, only_main_chain=False, main_chain='A'):
    chain_seq = []
    chain_dict = {}
    with open("/pubhome/hzhu02/GPSF/dataset/pdbbind_v2020/general_refine/"+pdb_id+"/"+pdb_id+"_protein.pdb", "r") as f:
        for line in f.readlines():
            if line.startswith("SEQRES"):
                line = line.strip().split()
                if line[2] not in chain_dict.keys():
                    chain_dict[line[2]]=''
                for i in range(4, len(line)):
                    aa = line[i]
                    if is_aa(aa, standard=True):
                        chain_dict[line[2]] += three_to_one(aa)
                    elif aa in {'HIE', 'HID'}:
                        chain_dict[line[2]] += 'H'
                    elif aa in {'CYX', 'CYM'}:
                        chain_dict[line[2]] += 'C'
                    else:
                        chain_dict[line[2]] += 'X' 
    if only_main_chain:
        try:
            chain_seq.append(chain_dict[main_chain])
        except:
            print(pdb_id)
    else:
        for key in chain_dict.keys():
            chain_seq.append(chain_dict[key])
    return chain_seq



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='create train/test sets for cross-validation separating by sequence similarity of protein targets and rdkit fingerprint similarity')
    parser.add_argument('--input',type=str,default="/pubhome/hzhu02/GPSF/dataset/pdbbind_v2020/index/pdbbind_v2020_refine.code",help="basic infomation of pdbid")
    parser.add_argument('--output_path',type=str, default="/pubhome/hzhu02/Redocked_pose/split_dataset/3_fold_ccv/general_refine",help='The prefix of the output file')  
    parser.add_argument('-c','--check',type=str,help='input name for folds to check for similarity')
    parser.add_argument('-n', '--number',type=int,default=3,help="number of folds to create/check. default=3")
    parser.add_argument('-s','--similarity',type=float,default=0.5,help='what percentage similarity to cluster by. default=0.5')
    parser.add_argument('-s2','--similarity_with_similar_ligand',type=float,default=0.4,help='what percentage similarity to cluster by when ligands are similar default=0.4')
    parser.add_argument('-l','--ligand_similarity',type=float,default=0.9,help='similarity threshold for ligands, default=0.9')
    parser.add_argument('--randomize',required=False,type=int,default=0,help='randomize inputs to get a different split, number is random seed')
    parser.add_argument('--main_chain_file', type=str, default="/pubhome/hzhu02/Redocked_pose/split_dataset/pfam/general/general_refine_main_chain_merge.csv")
    parser.add_argument('--only_main_chain', type=bool, default=True)
    args = parser.parse_args().__dict__

    threshold = 1 - args["similarity"] #similarity and distance are complementary
    threshold2 = 1 - args["similarity_with_similar_ligand"]
    ligand_threshold = args["ligand_similarity"] #this actually is a sim"]

    # calulate protein pairwise similarity
    # protein_seq = []
    # if args['only_main_chain']:
    #     main_chain_file = pd.read_csv(args['main_chain_file'],sep=",")
    #     pdbid_list = main_chain_file['pdb'].tolist()
    #     main_chain_list = main_chain_file['chain'].tolist()
    
        
        # for i in range(len(pdbid_list)):
        #     pdbid = pdbid_list[i]
        #     main_chain = main_chain_list[i]
        #     print(pdbid)
        #     seqs = get_chain(pdb_id=pdbid, only_main_chain=True, main_chain=main_chain)
        #     protein_seq.append(seqs)

    # else:
    #     complex_file=pd.read_csv(args["input"], header=None)
    #     complex_file.columns=['pdb_code']
    #     ##make sure 2e43 would not read as 2e00000000+6
    #     complex_file['pdb_code']=complex_file['pdb_code'].apply(lambda x: str(x)[0]+str(x).split("+")[0][-1]+str(x).split("+")[1] if len(str(x))>4 else str(x))
    
    #     ## get protein sequence similarity matrix

    #     pdbid_list = complex_file['pdb_code'].tolist()
    #     for pdbid in pdbid_list:
    #         seqs = get_chain(pdbid=pdbid)
    #         protein_seq.append(seqs)

    
    # print('Number of targets: {}'.format(len(pdbid_list)))


    # protein_smi_matrix = calProteinDistanceMatrix(protein_seq)
    # mol_fp_smi_matrix = calLigandSimilarity(fingerprints, pdbid_list)

    # pd.DataFrame(protein_smi_matrix, columns=pdbid_list).to_csv(args["output_path"] + "/protein_smi_only_main_chain_general_refine_pdbbind_2020.csv")

    ## calculate ligand fingerprint similarity




    # calculate ligand fingerprint similarity
    # fingerprint={}
    # for pdb in pdbid_list:
    #     print(pdb)
    #     # supplier = Chem.SDMolSupplier("/pubhome/hzhu02/GPSF/dataset/pdbbind_v2020/general_refine/"+pdb+"/"+pdb+"_ligand.sdf", sanitize=False, removeHs=False)
    #     # mol = supplier[0]
    #     mol = Chem.MolFromPDBFile("/pubhome/hzhu02/GPSF/dataset/pdbbind_v2020/general_refine/"+pdb+"/"+pdb+"_ligand_ob.pdb", sanitize=False, removeHs=False)
    #     fp = FingerprintMols.FingerprintMol(mol)
    #     fingerprint[pdb] = fp
    
    # ligand_similarity = calLigandSimilarity(fingerprint, pdbid_list)
    # pd.DataFrame(ligand_similarity, columns=pdbid_list).to_csv(args["output_path"] + "/ligand_similarity_pdbbind_general_refine_2020.csv", index=False)
    
    protein_dist = pd.read_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/3_fold_ccv/protein_smi_pdbbind_2020.csv")
    protein_dist.drop(columns=['Unnamed: 0'], inplace=True)
    pdb_list = protein_dist.columns.tolist()
    protein_smi_matrix = np.array(protein_dist)

    ligand_similarity = pd.read_csv("/pubhome/hzhu02/Redocked_pose/split_dataset/3_fold_ccv/ligand_similarity_pdbbind_2020.csv")
    ligand_similarity = np.array(ligand_similarity)
    if args["randomize"] != None:
        np.random.seed(args["randomize"])
        n = len(pdb_list)
        indices = np.random.choice(n,n,False)
        target_names = list(np.array(pdb_list)[indices])
        distanceMatrix = protein_smi_matrix[indices]
        distanceMatrix = distanceMatrix[:,indices]
        ligandsim = ligand_similarity[indices]
        ligandsim = ligandsim[:,indices]

    cluster_groups = calcClusterGroups(distanceMatrix, ligandsim, target_names, threshold, threshold2, ligand_threshold)
    print('{} clusters created'.format(len(cluster_groups)))
    
    for i, g in enumerate(cluster_groups):
        print('Cluster {}: {}'.format(i, ' '.join(str(t) for t in g)))

    print("Max cluster size: %d" % np.max([len(c) for c in cluster_groups]))

    # folds, foldmap = createFolds(cluster_groups, args["number"], args["randomize"])
    # for i, fold in enumerate(folds):
    #     print('{} targets in fold {}'.format(len(fold), i+1))
    #     write_file('%s%s'%(args["output_path"]+"/subset_", i+1), '\n'.join(fold))

    # folds, foldmap = loadFolds(args["check"], target_names, args["number"])
    # print('Checking {} train/test folds for {:.3f}% similarity'.format(args["check"], 100*args["similarity"]))
    # checkFolds(distanceMatrix, target_names, threshold, foldmap)








    
