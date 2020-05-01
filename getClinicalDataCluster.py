###################### get clinical data for each cluster #############################
# Author: T.J.M. Kuijpers data: 29-April-2020

# Input types to function
# input data: pandas dataframe containing the clinical data
# sampleColumn: the name of the column that contains the sample ids in the input data
# clusterData: 1d-array containing the sample ids
# columns to take: can be one string or 1d array containing multiple strings

def getclinicaldataclusters(inputData=None,sampleColumn=None,clusterData=None,columnsToTake=None):
    metaData=inputData[inputData[sampleColumn].isin(clusterData)]
    if columnsToTake==None:
        metaData = metaData
    else:
        metaData=metaData[:,columnsToTake]
    return metaData
