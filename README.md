# Fast Flow-based Random Walk with Restart in a Multi-query Setting
Codes of Fast Flow-based Random Walk with Restart in a Multi-query Setting

# Data
Data folder stores one dataset used in the paper, for the other datasets, download from https://snap.stanford.edu/data/

Data folder also stores the partition file of the corresponding graph datasets. Partition is obtained via Louvain method, which can be downloaded here: https://sites.google.com/site/findcommunities/

# Usage 
To use FlowR-ov and FlowR_HYB, cd to the corresponding folder and run the following respectively:

nohup matlab -nodisplay <runtime_overlap.m> log &

nohup matlab -nodisplay <runtime.m> log &

To change the datasets, change the names for the file and partition file in the runtime.m/runtime_overlap.m file, remember to change the directionality based on whether your graph is directed or undirected. 
# Paper:
For more details, pleae refer to our paper: http://web.eecs.umich.edu/~dkoutra/papers/18_FlowR_SDM_CR.pdf

For supplementary proofs, you can find it in the repository

To cite our paper:

@inproceedings{yan2018fast,
  title={Fast flow-based random walk with restart in a multi-query setting},
  
  author={Yan, Yujun and Heimann, Mark and Jin, Di and Koutra, Danai},
  
  booktitle={Proceedings of the 2018 SIAM International Conference on Data Mining},
  
  pages={342--350},
  
  year={2018},
  
  organization={SIAM}

