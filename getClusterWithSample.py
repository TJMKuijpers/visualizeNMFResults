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
    number=0
    for x in input:
        number=number+1
        if columnNames==None:
            clusterNumber="Cluster "+str(number)
            dataFrame=x.assign(Cluster=clusterNumber)
            dataFrameReturn.append(dataFrame)
        else:
             x=x.loc[:,columnNames]
             clusterNumber="Cluster "+str(number)
             dataFrame=x.assign(Cluster=clusterNumber)
             dataFrameReturn.append(dataFrame)
    dataFrameCluster=pd.concat(dataFrameReturn,axis=0)
    return dataFrameCluster

def meltdataframetolongformat(input=None,idVars=None,valueVars=None,varName=None):
    meltedFormat=pd.melt(input,id_vars=idVars,value_vars=valueVars,var_name=varName)
    return meltedFormat