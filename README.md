# Assessment of the Generalization Abilities of Machine-Learning Scoring Functions for Structure-Based Virtual Screening
Author: Hui Zhu, Jincai Yang*, and Niu Huang*
![](https://github.com/hnlab/generalization_benchmark/blob/main/png/Cover_Picture.png)

## Explanation
**split_dataset/cluster** collected scripts of Pocket Pfam-based clustering (Pfam-cluster) and Protein sequence-based clustering (Seq-cluster). *pdbbind_2020_cluster_result.csv* contained results of two clustering approaches


**split_dataset/3_fold** contained the training, validation and testing dataset for generalization ability benchmark in the paper.

**models/Descriptor_based_model**  contained source code of LR::V, LR::VR1, RF-Score, XGB::VR1 and NNScore. Other evaluated models were downloaded from individual paper.

|Models|Availability|
|--|--|
|Pafnucy|http://gitlab.com/cheminfIBB/pafnucy|
|OnionNet|http://github.com/zhenglz/onionnet/|
|SG-CNN|https://github.com/llnl/fast|
|IGN|https://github.com/zjujdj/InteractionGraphNet/tree/master|
|SIGN|https://github.com/PaddlePaddle/PaddleHelix/tree/dev/apps/drug_target_interaction/sign|
|GraphBAR|https://github.com/jtson82/graphbar|

**models/shap** is the Shapley Additive exPlanations (SHAP) analysis on RF-Score



## Citation
[Assessment of the Generalization Abilities of Machine-Learning Scoring Functions for Structure-Based Virtual Screening](https://pubs.acs.org/doi/10.1021/acs.jcim.2c01149)

