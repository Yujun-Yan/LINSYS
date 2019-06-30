# Fast Flow-based Random Walk with Restart in a Multi-query Setting
Codes of Fast Flow-based Random Walk with Restart in a Multi-query Setting
# Data folder stores one dataset used in the paper, for the other datasets, download from https://snap.stanford.edu/data/
# Data folder also stores the partition file of the corresponding graph datasets. Partition is obtained via Louvain method, which can be downloaded here: https://sites.google.com/site/findcommunities/
# To use FlowR-ov and FlowR_HYB, cd to the corresponding folder and run:
nohup matlab -nodisplay runtime_overlap.m > log &
nohup matlab -nodisplay runtime.m > log &
#respectively
#To change the datasets, change the names for the file and partition file in the runtime.m/runtime_overlap.m file, remember to change the directionality based on whether your graph is directed or undirected. 

