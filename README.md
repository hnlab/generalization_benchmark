# [Assessing the generalization abilities of machine-learning scoring functions for structure-based virtual screening](https://chemrxiv.org/engage/chemrxiv/article-details/62d53f80fe12e38913a7e287)
Author: Hui Zhu, Jincai Yang*, and Niu Huang*
![](https://github.com/hnlab/Generalization_benckmark/blob/main/plot_scripts/png/Cover_Picture.png)

## Explanation
**split_dataset/cluster** collected scripts of Pocket Pfam-based clustering (Pfam-cluster) and Protein sequence-based clustering (Seq-cluster). *pdbbind_2020_cluster_result.csv* contained results of two clustering approaches


**split_dataset/3_fold** contained the training, validation and testing dataset for generalization ability benchmark in the paper.

**models/Descriptor_based_model**  contained source code of LR::V, LR::VR1, RF-Score, XGB::VR1 and NNScore. Other evaluated models were downloaded from individual paper.

**models/shap** is the Shapley Additive exPlanations (SHAP) analysis on RF-Score



## Citation


