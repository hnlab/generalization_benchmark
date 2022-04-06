#!/bin/bash
#$ -N pair_smi
#$ -q ampere
##$ -l hostname=k231.hn.org
##$ -l gpu=1
#$ -pe ampere 18
#$ -o /pubhome/hzhu02/Redocked_pose/split_dataset/3_fold_ccv/general_refine/out
#$ -e /pubhome/hzhu02/Redocked_pose/split_dataset/3_fold_ccv/general_refine/error
#$ -wd /pubhome/hzhu02/Redocked_pose/split_dataset/3_fold_ccv/general_refine

conda activate fast
python /pubhome/hzhu02/Redocked_pose/split_dataset/3_fold_ccv/cluster.py