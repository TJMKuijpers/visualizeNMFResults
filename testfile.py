# Combine the two files to form one clinical data file
import pandas as pd
import numpy as np
import seaborn as sns
from getClusterWithSample import getclusterwithsample
from getExposureProfile import getexposureprofile
from createProfilePlot import createprofileplot
pathToCluster="C:/Users/tim.kuijpers/Desktop/EGM/Resultaten cohort totaal/cluster3.csv"
seperator=","
# Set the path to the data file that contains POP exposure profiles
pathToPOP="C:/Users/tim.kuijpers/Desktop/EGM/POPExposureObject_EGM.txt"
seperatorPopFile="\t"
# Set the path to the data file that contains the meta data (clinical info) of each sample
pathToMeta="C:/Users/tim.kuijpers/Desktop/EGM/ClinicalDataEGMSamples.csv"
seperatorMetaData=","
# Read the data
clusterData=pd.read_csv(pathToCluster,sep=seperator)
popExposuredata=pd.read_csv(pathToPOP,sep=seperatorPopFile)
clinicalDataSamples=pd.read_csv(pathToMeta,sep=seperatorMetaData)
# Get the shape of the data
clusterWithSamplesEGM = getclusterwithsample(clusterData,columnMembers='team',columnCluster='cluster')

popExposureCluster1=getexposureprofile(popExposuredata,clusterData=clusterWithSamplesEGM['cluster1'].team,idSampleColumn='egm_id')
createprofileplot(inputData=popExposureCluster1,title="Pop exposure cluster 1",xlabel="POPs",ylabel="Exposure[-]")