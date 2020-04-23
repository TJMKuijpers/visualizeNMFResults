import numpy as np
def getclusterwithsample(input = None,columnMembers = None,columnCluster = None):
    maxClusters = np.max(input[columnCluster])
    print(maxClusters)
    clustersInData=range(1,maxClusters+1)
    clustersWithSamples=dict()
    for x in clustersInData:
        samples=input[input[columnCluster]==x]
        clustersWithSamples["cluster"+str(x)]=samples
    return clustersWithSamples

def createdataframeclusters(input = None, columnMembers= None, columnCluster = None):
    return None