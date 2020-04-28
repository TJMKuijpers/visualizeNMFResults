import numpy as np
import pandas as pd
def getclusterwithsample(input = None,columnMembers = None,columnCluster = None):
    maxClusters = np.max(input[columnCluster])
    clustersInData=range(1,maxClusters+1)
    clustersWithSamples=dict()
    for x in clustersInData:
        samples=input[input[columnCluster]==x]
        clustersWithSamples["cluster"+str(x)]=samples
    return clustersWithSamples

def createdataframeclusters(input = None,columnNames=None):
    # First for each cluster create an extra column with "Cluster x "  as info
    dataFrameReturn=[]
    for y in range(1,input.__len__()+1):
        for x in input:
            if columnNames==None:
                clusterNumber="Cluster "+str(y)
                dataFrame=x.assign(Cluster=clusterNumber)
                dataFrameReturn.append(dataFrame)
            else:
                x=x.loc[:,columnNames]
                clusterNumber="Cluster "+str(y)
                dataFrame=x.assign(Cluster=clusterNumber)
                dataFrameReturn.append(dataFrame)
    dataFrameCluster=pd.concat(dataFrameReturn,axis=0)
    return dataFrameCluster
