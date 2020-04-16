# Function script to prepare names for radial tree

import pandas as pd

# Read the file with the cluster info
clusterInfo = pd.read_csv("E:/CCLE manusscript/Data/Results/Results used for paper/NMF clustering results/ClusterInfo_CCLE_8clusters.csv",sep=',',index_col=0)
# Read the meta data to get the histology classification for each sample
metaData=pd.read_csv("E:/CCLE manusscript/Data/Results/Results used for paper/Microsatellite_Instability_Info.csv",sep=";")
print(metaData.iloc[1:5,0:5])
# set the names for the radial tree
headParent="CCLE"
headChild="Cluster"
clusters=[1,2,3,4,5,6,7,8]

# Create the headParent nodes for the cell lines
centralParent=headParent+headChild
headParents=[headChild+' '+str(x) for x in clusters]

# Match the cell line id to the cell lineage
matchedCellLines=[metaData.CCLE_ID['Lineage'][metaData.CCLE_ID.str.contains(x)] for x in clusterInfo.team.iloc[x]]
print(matchedCellLines[1:4])
5637_URINARY_TRACT
