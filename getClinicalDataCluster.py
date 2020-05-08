###################### get clinical data for each cluster #############################
# Author: T.J.M. Kuijpers data: 29-April-2020

# Input types to function
# input data: pandas dataframe containing the clinical data
# sampleColumn: the name of the column that contains the sample ids in the input data
# clusterData: 1d-array containing the sample ids
# columns to take: can be one string or 1d array containing multiple strings

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def getclinicaldataclusters(inputData=None,sampleColumn=None,clusterData=None,columnsToTake=None):
    metaData=inputData[inputData[sampleColumn].isin(clusterData)]
    if columnsToTake==None:
        metaData = metaData
    else:
        metaData=metaData[:,columnsToTake]
    return metaData

def meltdataframetolongformat(input=None,idVars=None,valueVars=None,varName=None):
    meltedFormat=pd.melt(input,id_vars=idVars,value_vars=valueVars,var_name=varName)
    return meltedFormat

def plotGroupedBoxPlot(input=None,xVar=None,yVar=None,hueVar=None,title=None,xlabel=None,ylabel=None):
    boxplotGroup = sns.boxplot(x=xVar, y=yVar, hue=hueVar, data=input)
    boxplotGroup.axes.set_title(title, fontsize=16)
    boxplotGroup.set_xlabel(xlabel, fontsize=14)
    boxplotGroup.set_ylabel(ylabel, fontsize=14)
    boxplotGroup.tick_params(labelsize=10)
    plt.show()
    return boxplotGroup


def createCountDataFrame(input=None):
    countDict=dict()
    cluster=1
    for x in input:
        countDict['Cluster '+str(cluster)]=x.value_counts()
        cluster=cluster+1
    countDataFrame=pd.DataFrame.from_dict(countDict)
    return countDataFrame


def plotGroupedBarplot(input=None,xVar=None,yVar=None,hueVar=None,title=None,xlabel=None,ylabel=None):
    barplotGroup = sns.catplot(data=input,x=xVar,y=yVar,hue=hueVar,kind='count')
    barplotGroup.axes.set_title(title,fontsize=16)
    barplotGroup.set_xlabel(xlabel,fontsize=14)
    barplotGroup.set_ylabel(ylabel,fontsize=14)
    barplotGroup.tick_params(labelsize=10)
    plt.show()
    return barplotGroup