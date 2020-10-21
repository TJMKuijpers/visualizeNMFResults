import pandas as pd
import pickle
import numpy as np
from scipy.stats import zscore
import matplotlib.pyplot as plt
from getClusterWithSample import getclusterwithsample,createdataframeclusters,meltdataframetolongformat,determineClusterMembers,getclustermembersfromarray
from getExposureProfile import getexposureprofile
from scipy.stats import ttest_ind
from createProfilePlot import createprofileplot,creategroupprofileplot,createboxplotwithpvaluesbetweengroups
from calculateSignificantDifference import qqplotfordistribution,testnormaldistribution,nonparamtericmannwhitneyutest
from calculateCombinedZScore import calculateCombinedZScore
from extractFeaturesClusters import extractfeatures,createdisplotofscoredfeatures
from createCohortPlots import createCohortPlots
pathToExposureObject="E:/Phd Thesis/Chapter 4 subsection EGM/EGM/Genes and CpG methylation EGM/Data For ADA/POPExposureObject_EGM.txt"
seperator="\t"

# Read the exposure object
popExposure=pd.read_csv(pathToExposureObject,sep=seperator)
popExposure.index=popExposure['egm_id']
popExposure=popExposure.loc[:,['DDT', 'DDE', 'PCB.118', 'PCB.153', 'PCB.138', 'PCB.156','PCB.180', 'PCB.170', 'BDE.47']]
# calculate the Zscore
popExposureZScore=popExposure.apply(zscore)
popExposureZScore['egm_id']=popExposure.index.values
# Read the heavy metal exposure object
pathToMetalExposure="E:/Phd Thesis/Chapter 4 subsection EGM/EGM/Genes and CpG methylation EGM/Metals-v3.1.txt"
separatorMetals="\t"
metalExposure=pd.read_csv(pathToMetalExposure,sep=separatorMetals)
metalExposure.index=metalExposure['egm_id']
metalExposure=metalExposure.loc[:,["cadmium","lead"]]
metalExposureZscore=metalExposure.apply(zscore)
metalExposureZscore['egm_id']=metalExposure.index.values
# Read the meta data
metaDataPath="E:/Phd Thesis/Chapter 4 subsection EGM/EGM/Genes and CpG methylation EGM/MetadataObject.txt"
metaData=pd.read_csv(metaDataPath,sep="\t")
########################################## Analysis of total cohort exposure ##########################################################
# Analyse the exposure profiles of the clusters in the total cohort
"""
pathToCluster="E:/Phd Thesis/Chapter 4 subsection EGM/EGM/Resultaten cohort Italy/Cluster2Italy.csv"
clusterData=pd.read_csv(pathToCluster,sep=",")
clusterData.columns=['team','cluster']
clusterWithSamplesEGMTotalCohort = getclusterwithsample(clusterData,columnMembers='team',columnCluster='cluster')
dioxinlineZscore=calculateCombinedZScore(popExposureZScore,columnsToTake=['PCB.118','PCB.156','PCB.170'],axisToTake=1,scoreName='Dioxin-like PCBs')
nondioxinlineZscore=calculateCombinedZScore(popExposureZScore,columnsToTake=[ 'PCB.153', 'PCB.138', 'PCB.180'],axisToTake=1,scoreName='Nondioxin-like PCBs')
metalCombinedZscore=calculateCombinedZScore(metalExposure,columnsToTake=['cadmium','lead'],axisToTake=1,scoreName='Metals')
popsCombinedZScore=calculateCombinedZScore(popExposureZScore,columnsToTake=['DDE', 'PCB.118', 'PCB.153', 'PCB.138', 'PCB.156','PCB.180', 'PCB.170'],axisToTake=1,scoreName='combined POPs')
popExposureCluster1=getexposureprofile(popExposureZScore,clusterData=clusterWithSamplesEGMTotalCohort['cluster1'].team,idSampleColumn='egm_id')
popExposureCluster2=getexposureprofile(popExposureZScore,clusterData=clusterWithSamplesEGMTotalCohort['cluster2'].team,idSampleColumn='egm_id')
popExposureCluster3=getexposureprofile(popExposureZScore,clusterData=clusterWithSamplesEGMTotalCohort['cluster3'].team,idSampleColumn='egm_id')

popExposureCluster1_dlpcbs=getexposureprofile(dioxinlineZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster1'].team,idSampleColumn='egm_id')
popExposureCluster2_dlpcbs=getexposureprofile(dioxinlineZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster2'].team,idSampleColumn='egm_id')
popExposureCluster3_dlpcbs=getexposureprofile(dioxinlineZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster3'].team,idSampleColumn='egm_id')
popExposureCluster1_ndlpcbs=getexposureprofile(nondioxinlineZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster1'].team,idSampleColumn='egm_id')
popExposureCluster2_ndlpcbs=getexposureprofile(nondioxinlineZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster2'].team,idSampleColumn='egm_id')
popExposureCluster3_ndlpcbs=getexposureprofile(nondioxinlineZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster3'].team,idSampleColumn='egm_id')
popExposureCluster1_metals=getexposureprofile(metalCombinedZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster1'].team,idSampleColumn='egm_id')
popExposureCluster2_metals=getexposureprofile(metalCombinedZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster2'].team,idSampleColumn='egm_id')
popExposureCluster3_metals=getexposureprofile(metalCombinedZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster3'].team,idSampleColumn='egm_id')
popExposureCluster1_pops=getexposureprofile(popsCombinedZScore,clusterData=clusterWithSamplesEGMTotalCohort['cluster1'].team,idSampleColumn='egm_id')
popExposureCluster2_pops=getexposureprofile(popsCombinedZScore,clusterData=clusterWithSamplesEGMTotalCohort['cluster2'].team,idSampleColumn='egm_id')
popExposureCluster3_pops=getexposureprofile(popsCombinedZScore,clusterData=clusterWithSamplesEGMTotalCohort['cluster3'].team,idSampleColumn='egm_id')

combinedDataFrame_dlpcbs=createdataframeclusters(input=[popExposureCluster1_dlpcbs,popExposureCluster2_dlpcbs,popExposureCluster3_dlpcbs],columnNames=['Zscore'])
#combinedDataFrame_ndlpcbs=createdataframeclusters(input=[popExposureCluster1_ndlpcbs,popExposureCluster2_ndlpcbs,popExposureCluster3_ndlpcbs],columnNames=['Zscore'])
#combinedDataFrame_metals=createdataframeclusters(input=[popExposureCluster1_metals,popExposureCluster2_metals,popExposureCluster3_metals],columnNames=['Zscore'])
#combinedDataFrame_pops=createdataframeclusters(input=[popExposureCluster1_pops,popExposureCluster2_pops,popExposureCluster3_pops],columnNames=['Zscore'])
#combinedDataFrameLong_ndlpcbs=meltdataframetolongformat(input=combinedDataFrame_dlpcbs,idVars=['Cluster'],valueVars=['Zscore'],varName='DL-PCBs')
#combinedDataFrameLong_dlpcbs=meltdataframetolongformat(input=combinedDataFrame_ndlpcbs,idVars=['Cluster'],valueVars=['Zscore'],varName='NDL-PCBs')
#combinedDataFrameLong_metals=meltdataframetolongformat(input=combinedDataFrame_metals,idVars=['Cluster'],valueVars=['Zscore'],varName='Metals')
#combinedDataFrameLong_pops=meltdataframetolongformat(input=combinedDataFrame_pops,idVars=['Cluster'],valueVars=['Zscore'],varName='POPs')
#dlpcbsPlot=creategroupprofileplot(inputData=combinedDataFrameLong_dlpcbs,xVar='Cluster',yVar='value',hueVar='Cluster',title='Total Cohort',xlabel='DL-PCBs',ylabel='Exposure [Z-score]')
#ndlpcbsPlot=creategroupprofileplot(inputData=combinedDataFrameLong_ndlpcbs,xVar='Cluster',yVar='value',hueVar='Cluster',title='Total Cohort',xlabel='NDL-PCBs',ylabel='Exposure [Z-score]')
#metalsPlot=creategroupprofileplot(inputData=combinedDataFrameLong_metals,xVar='Cluster',yVar='value',hueVar='Cluster',title='Total Cohort',xlabel='Metals',ylabel='Exposure [Z-score]')
#popsPlot=creategroupprofileplot(inputData=combinedDataFrameLong_pops,xVar='Cluster',yVar='value',hueVar='Cluster',title='Total Cohort',xlabel='POPs',ylabel='Exposure [Z-score]')
# Create a grouped boxplot of the exposures
#combinedDataFramePops_DDEDDT=createdataframeclusters(input=[popExposureCluster1,popExposureCluster2,popExposureCluster3],columnNames=['DDT', 'DDE'])
#combinedDataFramePopLong_DDEDDT=meltdataframetolongformat(input=combinedDataFramePops_DDEDDT,idVars=['Cluster'],valueVars=['DDT', 'DDE'],varName='POP')
#plotDDEDDTProfileGroup=creategroupprofileplot(inputData=combinedDataFramePopLong_DDEDDT,xVar='POP',yVar='value',hueVar='Cluster',title='Total Cohort',xlabel='POPs',ylabel='Exposure [Z-score]')
#combinedDataFramePops_PCBs=createdataframeclusters(input=[popExposureCluster1,popExposureCluster2,popExposureCluster3],columnNames=['PCB.118', 'PCB.153', 'PCB.138', 'PCB.156','PCB.180', 'PCB.170'])
#combinedDataFramePopLong_PCBs=meltdataframetolongformat(input=combinedDataFramePops_PCBs,idVars=['Cluster'],valueVars=['PCB.118', 'PCB.153', 'PCB.138', 'PCB.156','PCB.180', 'PCB.170'],varName='POP')
#plotPCBsProfileGroup=creategroupprofileplot(inputData=combinedDataFramePopLong_PCBs,xVar='POP',yVar='value',hueVar='Cluster',title='Total Cohort',xlabel='POPs',ylabel='Exposure [Z-score]')
#Group1VersusGroup2DDETotal=nonparamtericmannwhitneyutest(popExposureCluster1.loc[:,'DDE'],popExposureCluster2.loc[:,'DDE'])
#Group1VersusGroup3DDETotal=nonparamtericmannwhitneyutest(popExposureCluster1.loc[:,'DDE'],popExposureCluster3.loc[:,'DDE'])
#Group1VersusGroup3DDETotal=nonparamtericmannwhitneyutest(popExposureCluster2.loc[:,'DDE'],popExposureCluster3.loc[:,'DDE'])
#print("Cluster 1 vs cluster 2 MWU test DDE Total",Group1VersusGroup2DDETotal)
#print("Cluster 1 vs cluster 3 MWU test DDE Total",Group1VersusGroup3DDETotal)
#print("Cluster 2 vs cluster 3 MWU test DDE Total",Group1VersusGroup3DDETotal)
#DDETotalCohortPlot=createboxplotwithpvaluesbetweengroups(combinedDataFramePopLong_DDEDDT,xVar='POP',yVar='value',hueVar='Cluster',title='',xlabel='Cluster',ylabel='Exposure [Z-score]',significanceColumns=[[0.75,1.25],[0.75,1],],columnForMax='value')
# Create the profile plots for the heavy metals
#metalExposureCluster1=getexposureprofile(metalExposureZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster1'].team,idSampleColumn='egm_id')
#metalExposureCluster2=getexposureprofile(metalExposureZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster2'].team,idSampleColumn='egm_id')
#metalExposureCluster3=getexposureprofile(metalExposureZscore,clusterData=clusterWithSamplesEGMTotalCohort['cluster3'].team,idSampleColumn='egm_id')
#metalExposureProfilesCombinedTotalCohort=createdataframeclusters(input=[metalExposureCluster1,metalExposureCluster2,metalExposureCluster3],columnNames=["cadmium","lead"])
#combinedDataFrameMetalsLong_TotalCohort=meltdataframetolongformat(input=metalExposureProfilesCombinedTotalCohort,idVars=['Cluster'],valueVars=["cadmium","lead"],varName='metals')
#plotMetalsTotalCohortProfileGroup=creategroupprofileplot(inputData=combinedDataFrameMetalsLong_TotalCohort,xVar='metals',yVar='value',hueVar='Cluster',title='Total Cohort',xlabel='Metals',ylabel='Exposure [Z-score]')

#Group1VersusGroup2LeadTotal=nonparamtericmannwhitneyutest(metalExposureCluster1.loc[:,'lead'],metalExposureCluster2.loc[:,'lead'])
#Group1VersusGroup3LeadTotal=nonparamtericmannwhitneyutest(metalExposureCluster1.loc[:,'lead'],metalExposureCluster3.loc[:,'lead'])
#Group2VersusGroup3LeadTotal=nonparamtericmannwhitneyutest(metalExposureCluster2.loc[:,'lead'],metalExposureCluster3.loc[:,'lead'])
#print("Cluster 1 vs cluster 2 MWU test Lead Total",Group1VersusGroup2LeadTotal)
#print("Cluster 1 vs cluster 3 MWU test  Lead Total",Group1VersusGroup3LeadTotal)
#print("Cluster 2 vs cluster 3 MWU Lead test Total",Group2VersusGroup3LeadTotal)
#leadTotalCohortPlot=createboxplotwithpvaluesbetweengroups(combinedDataFrameMetalsLong_TotalCohort,xVar='metals',yVar='value',hueVar='Cluster',title='',xlabel='Cluster',ylabel='Exposure [Z-score]',significanceColumns=[[0.75,1.25],[0.75,1],[1,1.25]],columnForMax='value')
#Group1VersusGroup2CadmiumTotal=nonparamtericmannwhitneyutest(metalExposureCluster1.loc[:,'cadmium'],metalExposureCluster2.loc[:,'cadmium'])
#Group1VersusGroup3CadmiumTotal=nonparamtericmannwhitneyutest(metalExposureCluster1.loc[:,'cadmium'],metalExposureCluster3.loc[:,'cadmium'])
#Group2VersusGroup3CadmiumTotal=nonparamtericmannwhitneyutest(metalExposureCluster2.loc[:,'cadmium'],metalExposureCluster3.loc[:,'cadmium'])
#print("Cluster 1 vs cluster 2 MWU test cadmium Total",Group1VersusGroup2CadmiumTotal)
#print("Cluster 1 vs cluster 3 MWU test cadmium Total",Group1VersusGroup3CadmiumTotal)
#print("Cluster 2 vs cluster 3 MWU test cadmium Total",Group2VersusGroup3CadmiumTotal)
"""
############################################# Create the profile plots for the Swedish cohort ###############################################################
# Extract the clusters from the H coefficient matrix

pathToHMatrix="E:/Phd Thesis/Chapter 4 subsection EGM/EGM/Resultaten cohort Sweden/H_k3_100runs.pkl"
pathToGeneExpression="E:/Phd Thesis/Chapter 4 subsection EGM/EGM/Genes and CpG methylation EGM/Data For ADA/Cohorts/GeneExpression_Sweden.csv"
geneExpressionData=pd.read_csv(pathToGeneExpression,sep=",",index_col=0)
loadPickleObject=open(pathToHMatrix,"rb")
matrixHSolution=pickle.load(loadPickleObject)
print("Number of simulations of H is %d" % matrixHSolution[0].__len__())
matrixH=matrixHSolution[0][9]
clusters,output=determineClusterMembers(matrixH, geneExpressionData.columns.values)
cluster1Sweden=clusters[0]
cluster2Sweden=clusters[1]
cluster3Sweden=clusters[2]
clustersSweden=getclustermembersfromarray(input=[cluster1Sweden,cluster2Sweden,cluster3Sweden])
clustersSweden.columns=['team','cluster']
# Calculate POP exposure Z score for only the swedish cohort
popExposureSweden=popExposure.loc[popExposure.index.isin(clustersSweden.team),:]
popExposureSwedenZScore=popExposureSweden.apply(zscore)
popExposureSwedenZScore['egm_id']=popExposureSweden.index.values
# calculate metal exposure zscore for the swedish cohort
metalExposureSweden=metalExposure.loc[metalExposure.index.isin(clustersSweden.team),:]
metalExposureSwedenZScore=metalExposureSweden.apply(zscore)
metalExposureSwedenZScore['egm_id']=metalExposureSweden.index.values

popExposureCluster1Sweden=getexposureprofile(popExposureSwedenZScore,clusterData=clustersSweden[clustersSweden['cluster']=='cluster 1'].team,idSampleColumn='egm_id')
popExposureCluster2Sweden=getexposureprofile(popExposureSwedenZScore,clusterData=clustersSweden[clustersSweden['cluster']=='cluster 2'].team,idSampleColumn='egm_id')
popExposureCluster3Sweden=getexposureprofile(popExposureSwedenZScore,clusterData=clustersSweden[clustersSweden['cluster']=='cluster 3'].team,idSampleColumn='egm_id')
combinedDataFramePops_DDEDDT_Sweden=createdataframeclusters(input=[popExposureCluster1Sweden,popExposureCluster1Sweden,popExposureCluster1Sweden],columnNames=['DDT', 'DDE'])
combinedDataFramePopLong_DDEDDT_Sweden=meltdataframetolongformat(input=combinedDataFramePops_DDEDDT_Sweden,idVars=['Cluster'],valueVars=['DDT', 'DDE'],varName='POP')
plotDDEDDTProfileGroupSweden=creategroupprofileplot(inputData=combinedDataFramePopLong_DDEDDT_Sweden,xVar='POP',yVar='value',hueVar='Cluster',title='Swedish cohort',xlabel='POPs',ylabel='Exposure [Z-score]')
combinedDataFramePops_PCBs_Sweden=createdataframeclusters(input=[popExposureCluster1Sweden,popExposureCluster2Sweden,popExposureCluster3Sweden],columnNames=['PCB.118', 'PCB.153', 'PCB.138', 'PCB.156','PCB.180', 'PCB.170'])
combinedDataFramePopLong_PCBs_Sweden=meltdataframetolongformat(input=combinedDataFramePops_PCBs_Sweden,idVars=['Cluster'],valueVars=['PCB.118', 'PCB.153', 'PCB.138', 'PCB.156','PCB.180', 'PCB.170'],varName='POP')
plotPCBsProfileGroup_Sweden=creategroupprofileplot(inputData=combinedDataFramePopLong_PCBs_Sweden,xVar='POP',yVar='value',hueVar='Cluster',title='Swedish cohort',xlabel='POPs',ylabel='Exposure [Z-score]')
# Test for significant difference
normalityDDEcluster1=testnormaldistribution(popExposureCluster1Sweden.loc[:,'DDE'])
print("Normality cluster 1: ", normalityDDEcluster1)
normalityDDEcluster2=testnormaldistribution(popExposureCluster2Sweden.loc[:,'DDE'])
print("Normality cluster 2: ", normalityDDEcluster2)
normalityDDEcluster3=testnormaldistribution(popExposureCluster3Sweden.loc[:,'DDE'])
print("Normality cluster 3: ", normalityDDEcluster3)
Group1VersusGroup2=nonparamtericmannwhitneyutest(popExposureCluster1Sweden.loc[:,'DDE'],popExposureCluster2Sweden.loc[:,'DDE'])
Group1VersusGroup3=nonparamtericmannwhitneyutest(popExposureCluster1Sweden.loc[:,'DDE'],popExposureCluster3Sweden.loc[:,'DDE'])
Group2VersusGroup3=nonparamtericmannwhitneyutest(popExposureCluster2Sweden.loc[:,'DDE'],popExposureCluster3Sweden.loc[:,'DDE'])
print("Cluster 1 vs cluster 2 MWU test DDE Sweden",Group1VersusGroup2)
print("Cluster 1 vs cluster 3 MWU test DDE Sweden",Group1VersusGroup3)
print("Cluster 2 vs cluster 3 MWU test DDE Sweden",Group2VersusGroup3)

# Metal exposure plot
metalExposureCluster1Sweden=getexposureprofile(metalExposureSwedenZScore,clusterData=clustersSweden[clustersSweden['cluster']=='cluster 1'].team,idSampleColumn='egm_id')
metalExposureCluster2Sweden=getexposureprofile(metalExposureSwedenZScore,clusterData=clustersSweden[clustersSweden['cluster']=='cluster 2'].team,idSampleColumn='egm_id')
metalExposureCluster3Sweden=getexposureprofile(metalExposureSwedenZScore,clusterData=clustersSweden[clustersSweden['cluster']=='cluster 3'].team,idSampleColumn='egm_id')
metalExposureProfilesCombinedSweden=createdataframeclusters(input=[metalExposureCluster1Sweden,metalExposureCluster2Sweden,metalExposureCluster3Sweden],columnNames=["cadmium","lead"])
combinedDataFrameMetalsLong_Sweden=meltdataframetolongformat(input=metalExposureProfilesCombinedSweden,idVars=['Cluster'],valueVars=["cadmium","lead"],varName='metals')
plotMetalsSwedenProfileGroup=creategroupprofileplot(inputData=combinedDataFrameMetalsLong_Sweden,xVar='metals',yVar='value',hueVar='Cluster',title='Swedish Cohort',xlabel='Metals',ylabel='Exposure [Z-score]')
Group1VersusGroup2Lead=nonparamtericmannwhitneyutest(metalExposureCluster1Sweden.loc[:,'lead'],metalExposureCluster2Sweden.loc[:,'lead'])
Group1VersusGroup3Lead=nonparamtericmannwhitneyutest(metalExposureCluster1Sweden.loc[:,'lead'],metalExposureCluster3Sweden.loc[:,'lead'])
Group2VersusGroup3Lead=nonparamtericmannwhitneyutest(metalExposureCluster2Sweden.loc[:,'lead'],metalExposureCluster3Sweden.loc[:,'lead'])
print("Cluster 1 vs cluster 2 MWU test Lead Sweden",Group1VersusGroup2Lead)
print("Cluster 1 vs cluster 3 MWU test Lead Sweden",Group1VersusGroup3Lead)
print("Cluster 2 vs cluster 3 MWU test Lead Sweden",Group2VersusGroup3Lead)

Group1VersusGroup2Cadmium=nonparamtericmannwhitneyutest(metalExposureCluster1Sweden.loc[:,'cadmium'],metalExposureCluster2Sweden.loc[:,'cadmium'])
Group1VersusGroup3Cadmium=nonparamtericmannwhitneyutest(metalExposureCluster1Sweden.loc[:,'cadmium'],metalExposureCluster3Sweden.loc[:,'cadmium'])
Group2VersusGroup3Cadmium=nonparamtericmannwhitneyutest(metalExposureCluster2Sweden.loc[:,'cadmium'],metalExposureCluster3Sweden.loc[:,'cadmium'])
print("Cluster 1 vs cluster 2 MWU test cadmium Sweden",Group1VersusGroup2Cadmium)
print("Cluster 1 vs cluster 3 MWU test cadmium Sweden",Group1VersusGroup3Cadmium)
print("Cluster 2 vs cluster 3 MWU test cadmium Sweden",Group2VersusGroup3Cadmium)

#leadSwedenPlot=createboxplotwithpvaluesbetweengroups(combinedDataFrameMetalsLong_Sweden,xVar='metals',yVar='value',hueVar='Cluster',title='',xlabel='Cluster',ylabel='Exposure [Z-score]',significanceColumns=[[0.75,1.25],[1,1.25]],columnForMax='value')
####################################################### Create the profile plots for the Italian cohort ##############################################
"""
pathToClusterItaly="E:/Phd Thesis/Chapter 4 subsection EGM/EGM/Resultaten cohort Italy/Cluster2Italy.csv"
separatorItaly=","
clusterDataItaly=pd.read_csv(pathToClusterItaly,sep=separatorItaly)
# Calculate the POP exposure zscore for italy
popExposureItaly=popExposure.loc[popExposure.index.isin(clusterDataItaly.team),:]
popExposureItalyZscore=popExposureItaly
popExposureItalyZscore=popExposureItaly.apply(zscore)
popExposureItalyZscore['egm_id']=popExposureItaly.index.values
# Calculate metal exposure z score for italy
metalExposureItaly=metalExposure.loc[metalExposure.index.isin(clusterDataItaly.team),:]
metalExposureItalyZScore=metalExposureItaly.apply(zscore)
metalExposureItalyZScore['egm_id']=metalExposureItaly.index.values

popExposureCluster1Italy=getexposureprofile(popExposureItalyZscore,clusterData=clusterDataItaly[clusterDataItaly['cluster']==1].team,idSampleColumn='egm_id')
popExposureCluster2Italy=getexposureprofile(popExposureItalyZscore,clusterData=clusterDataItaly[clusterDataItaly['cluster']==2].team,idSampleColumn='egm_id')
#popExposureCluster3Italy=getexposureprofile(popExposureItalyZscore,clusterData=clusterDataItaly[clusterDataItaly['cluster']==3].team,idSampleColumn='egm_id')
combinedDataFramePops_DDEDDT_Italy=createdataframeclusters(input=[popExposureCluster1Italy,popExposureCluster2Italy],columnNames=['DDT', 'DDE'])
combinedDataFramePopLong_DDEDDT_Italy=meltdataframetolongformat(input=combinedDataFramePops_DDEDDT_Italy,idVars=['Cluster'],valueVars=['DDE'],varName='POP')
plotDDEDDTProfileGroupItaly=creategroupprofileplot(inputData=combinedDataFramePopLong_DDEDDT_Italy,xVar='POP',yVar='value',hueVar='Cluster',title='Combined Cohorts',xlabel='POPs',ylabel='Exposure [Z-score]')
combinedDataFramePops_PCBs_Italy=createdataframeclusters(input=[popExposureCluster1Italy,popExposureCluster2Italy],columnNames=['PCB.118', 'PCB.153', 'PCB.138', 'PCB.156','PCB.180', 'PCB.170'])
combinedDataFramePopLong_PCBs_Italy=meltdataframetolongformat(input=combinedDataFramePops_PCBs_Italy,idVars=['Cluster'],valueVars=['PCB.118', 'PCB.153', 'PCB.138', 'PCB.156','PCB.180', 'PCB.170'],varName='POP')
plotPCBsProfileGroup_Italy=creategroupprofileplot(inputData=combinedDataFramePopLong_PCBs_Italy,xVar='POP',yVar='value',hueVar='Cluster',title='Combined Cohorts',xlabel='POPs',ylabel='Exposure [Z-score]')
Group1VersusGroup2=nonparamtericmannwhitneyutest(popExposureCluster1Italy.loc[:,'DDE'],popExposureCluster1Italy.loc[:,'DDE'])
#print("Cluster 1 vs cluster 2 MWU test DDE Italy",Group1VersusGroup2)
dioxinlineZscoreItaly=calculateCombinedZScore(popExposureItalyZscore,columnsToTake=['PCB.118','PCB.156','PCB.170'],axisToTake=1,scoreName='Dioxin-like PCBs')
nondioxinlineZscoreItaly=calculateCombinedZScore(popExposureItalyZscore,columnsToTake=[ 'PCB.153', 'PCB.138', 'PCB.180'],axisToTake=1,scoreName='Nondioxin-like PCBs')
popExposureCluster1Italy_dlpcbs=getexposureprofile(dioxinlineZscoreItaly,clusterData=clusterDataItaly[clusterDataItaly['cluster']==1].team,idSampleColumn='egm_id')
popExposureCluster2Italy_dlpcbs=getexposureprofile(dioxinlineZscoreItaly,clusterData=clusterDataItaly[clusterDataItaly['cluster']==2].team,idSampleColumn='egm_id')

popExposureCluster1Italy_ndlpcbs=getexposureprofile(nondioxinlineZscoreItaly,clusterData=clusterDataItaly[clusterDataItaly['cluster']==1].team,idSampleColumn='egm_id')
popExposureCluster2Italy_ndlpcbs=getexposureprofile(nondioxinlineZscoreItaly,clusterData=clusterDataItaly[clusterDataItaly['cluster']==2].team,idSampleColumn='egm_id')
combinedDataFrameItaly_dlpcbs=createdataframeclusters(input=[popExposureCluster1Italy_dlpcbs,popExposureCluster2Italy_dlpcbs],columnNames=['Zscore'])
combinedDataFrameItaly_ndlpcbs=createdataframeclusters(input=[popExposureCluster1Italy_ndlpcbs,popExposureCluster2Italy_ndlpcbs],columnNames=['Zscore'])

combinedDataFrameLongItaly_ndlpcbs=meltdataframetolongformat(input=combinedDataFrameItaly_dlpcbs,idVars=['Cluster'],valueVars=['Zscore'],varName='DL-PCBs')
combinedDataFrameLongItaly_dlpcbs=meltdataframetolongformat(input=combinedDataFrameItaly_ndlpcbs,idVars=['Cluster'],valueVars=['Zscore'],varName='NDL-PCBs')
dlpcbsPlot=creategroupprofileplot(inputData=combinedDataFrameLongItaly_ndlpcbs,xVar='Cluster',yVar='value',hueVar='Cluster',title='Combined Cohorts',xlabel='DL-PCBs',ylabel='Exposure [Z-score]')
ndlpcbsPlot=creategroupprofileplot(inputData=combinedDataFrameLongItaly_dlpcbs,xVar='Cluster',yVar='value',hueVar='Cluster',title='Combined Cohorts',xlabel='NDL-PCBs',ylabel='Exposure [Z-score]')

popsCombinedZScore=calculateCombinedZScore(popExposureZScore,columnsToTake=['DDE', 'PCB.118', 'PCB.153', 'PCB.138', 'PCB.156','PCB.180', 'PCB.170'],axisToTake=1,scoreName='combined POPs')
popExposureCluster1Italy_pops=getexposureprofile(popsCombinedZScore,clusterData=clusterDataItaly[clusterDataItaly['cluster']==1].team,idSampleColumn='egm_id')
popExposureCluster2Italy_pops=getexposureprofile(popsCombinedZScore,clusterData=clusterDataItaly[clusterDataItaly['cluster']==2].team,idSampleColumn='egm_id')
combinedDataFrame_pops=createdataframeclusters(input=[popExposureCluster1Italy_pops,popExposureCluster2Italy_pops],columnNames=['Zscore'])
combinedDataFrameLong_pops=meltdataframetolongformat(input=combinedDataFrame_pops,idVars=['Cluster'],valueVars=['Zscore'],varName='POPs')
popsPlot=creategroupprofileplot(inputData=combinedDataFrameLong_pops,xVar='Cluster',yVar='value',hueVar='Cluster',title='Combined Cohorts',xlabel='POPs',ylabel='Exposure [Z-score]')

# Metal exposure plot
metalExposureCluster1Italy=getexposureprofile(metalExposureItalyZScore,clusterData=clusterDataItaly[clusterDataItaly['cluster']==1].team,idSampleColumn='egm_id')
metalExposureCluster2Italy=getexposureprofile(metalExposureItalyZScore,clusterData=clusterDataItaly[clusterDataItaly['cluster']==2].team,idSampleColumn='egm_id')
#metalExposureCluster3Italy=getexposureprofile(metalExposureItalyZScore,clusterData=clusterDataItaly[clusterDataItaly['cluster']==3].team,idSampleColumn='egm_id')
metalExposureProfilesCombinedItaly=createdataframeclusters(input=[metalExposureCluster1Italy,metalExposureCluster2Italy],columnNames=["cadmium","lead"])
combinedDataFrameMetalsLong_Italy=meltdataframetolongformat(input=metalExposureProfilesCombinedItaly,idVars=['Cluster'],valueVars=["cadmium","lead"],varName='metals')
plotMetalsItalyProfileGroup=creategroupprofileplot(inputData=combinedDataFrameMetalsLong_Italy,xVar='metals',yVar='value',hueVar='Cluster',title='Combined Cohorts',xlabel='Metals',ylabel='Exposure [Z-score]')

Group1VersusGroup2=nonparamtericmannwhitneyutest(metalExposureCluster1Italy.loc[:,'lead'],metalExposureCluster2Italy.loc[:,'lead'])
print("Cluster 1 vs cluster 2 MWU test Lead Italy",Group1VersusGroup2)
Group1VersusGroup2=nonparamtericmannwhitneyutest(metalExposureCluster1Italy.loc[:,'cadmium'],metalExposureCluster2Italy.loc[:,'cadmium'])
print("Cluster 1 vs cluster 2 MWU test cadmium Italy",Group1VersusGroup2)
"""
#metalPlotIlatySignificant=createboxplotwithpvaluesbetweengroups(combinedDataFrameMetalsLong_Italy,xVar='metals',yVar='value',hueVar='Cluster',title='',xlabel='Cluster',ylabel='Exposure [Z-score]',significanceColumns=[[0.75,1.25]],columnForMax='value')
# Note popExposureCluster2Italy is still the old name for total cohort
"""
HmatrixPath="E:/Project Envirogenomarkers/Sex Corrected Results/TotalCohortSimulation_SexCorrected_Hmatrix.csv"
WmatrixGenePath="E:/Project Envirogenomarkers/Sex Corrected Results/TotalCohortSimulation_SexCorrected_WmatrixGene.csv"
WMatrixCpGPath="E:/Project Envirogenomarkers/Sex Corrected Results/TotalCohortSimulation_SexCorrected_WmatrixCpG.csv"
Hmatrix=pd.read_csv(HmatrixPath,sep=",",index_col=0)
WmatrixGene=pd.read_csv(WmatrixGenePath,sep=",",index_col=0)
WmatrixCpG=pd.read_csv(WMatrixCpGPath,sep=",",index_col=0)
# Get the clustering from the H matrix
cluster1=[]
cluster2=[]
for x in range(Hmatrix.shape[1]):
    minpos=np.argmin(Hmatrix.iloc[:,x].values)
    if minpos==0:
        cluster1.append(Hmatrix.columns[x])
    if minpos==1:
        cluster2.append(Hmatrix.columns[x])
Cluster1=pd.DataFrame({'Sample':cluster1})
Cluster2=pd.DataFrame({'Sample':cluster2})

scoredMatrixCpGs,selCpGs,featuresScoredCpGs=extractfeatures(WmatrixCpG.to_numpy(),methodScore='Kim',score_threshold=None,clustersize=2)
scoredMatrixGenes,selGenes,featuresScoredGenes=extractfeatures(WmatrixGene.to_numpy(),methodScore='Kim',score_threshold=None,clustersize=2)
scoredMatrixCpGsMax,selCpGsMax,featuresScoredCpGsMax=extractfeatures(WmatrixCpG.to_numpy(),methodScore='Max_Kim_Score',score_threshold=0.7,clustersize=2)
scoredMatrixGenesMax,selGenesMax,featuresScoredGenesMax=extractfeatures(WmatrixGene.to_numpy(),methodScore='Max_Kim_Score',score_threshold=0.7,clustersize=2)

pathToGeneExpression="E:/Phd Thesis/Chapter 4 subsection EGM/EGM/Genes and CpG methylation EGM/Data For ADA/GeneExpressionTotalCohorts.txt"
pathToMethylationData="E:/Phd Thesis/Chapter 4 subsection EGM/EGM/Genes and CpG methylation EGM/Data For ADA/methylationDataObject.txt"
geneExpressionData=pd.read_csv(pathToGeneExpression,sep="\t",index_col=0)
methylationData=pd.read_csv(pathToMethylationData,sep="\t")
geneProbesCluster1=WmatrixGene.index[featuresScoredGenes[0]]
geneProbesCluster2=WmatrixGene.index[featuresScoredGenes[1]]
cpgVariablesCluster1=WmatrixCpG.index[featuresScoredCpGsMax[0]]
cpgVariablesCluster2=WmatrixCpG.index[featuresScoredCpGsMax[1]]
import seaborn as sns
testProbes=geneProbesCluster1[1:15]
geneExpressionCluster1=geneExpressionData.loc[geneExpressionData.index.isin(testProbes),geneExpressionData.columns.isin(Cluster1.Sample)]
geneExpressionCluster1=geneExpressionCluster1.T.assign(Cluster='Cluster 1')
geneExpressionCluster2=geneExpressionData.loc[geneExpressionData.index.isin(testProbes),geneExpressionData.columns.isin(Cluster2.Sample)]
geneExpressionCluster2=geneExpressionCluster2.T.assign(Cluster='Cluster 2')
genetest=pd.concat([geneExpressionCluster1,geneExpressionCluster2])
genemelt=pd.melt(genetest,id_vars=['Cluster'], var_name=['Probe'])
plt.figure()
plt.title('gene expression')
sns.boxplot(x='Probe',y='value',hue='Cluster',data=genemelt)
plt.xticks(rotation=90)
plt.show()

"""
