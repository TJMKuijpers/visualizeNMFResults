import numpy as np
import pandas as pd

def getconsensusmatrixfromhmatrix(input=None,sampleNames=None):
    consensusmatrixH=consensusmatrix(Hdata=input, samplenames=sampleNames)
    return consensusmatrixH

def consensusMatrix(Hdata, samplenames):
    CoAssociationMatrices = []
    for x in range(len(Hdata)):
        MembersOfClusters, occurenceCluster = determineClusterMembers(Hdata[x], samplenames)
        CoAssociationMatrix = coAssociationMatrix(occurenceCluster)
        CoAssociationMatrices.append(CoAssociationMatrix)
        Consensusmatrix = sum(CoAssociationMatrices)
        Consensusmatrix = Consensusmatrix / len(Hdata)

    return Consensusmatrix


def coAssociationMatrix(clusterlist):
    results = [[int(x == y) for y in clusterlist] for x in clusterlist]
    coAssociation = np.array(results)
    return coAssociation

def determineClusterMembers(Hmatrix, samplenames):
    "For each column in H, maximum value will be determined and right cluster assigned"
    index_maxvalue = np.argmax(Hmatrix, axis=0)
    ClusterMembers = []
    for cluster in range(np.min(index_maxvalue), np.max(index_maxvalue) + 1):
        ClusterMembers.append([i for indx, i in enumerate(samplenames) if index_maxvalue[indx] == cluster])
    return ClusterMembers, index_maxvalue

def getclustermembersfromarray(input=None):
    dataframeClusters=[]
    dataFrameList=[]
    k=1
    for x in input:
        tmpDataFrame=pd.DataFrame({'Sample':x})
        tmpDataFrame=tmpDataFrame.assign(Cluster='cluster '+str(k))
        dataFrameList.append(tmpDataFrame)
        k=k+1
    dataframeClusters=pd.concat(dataFrameList,axis=0)
    return dataframeClusters


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