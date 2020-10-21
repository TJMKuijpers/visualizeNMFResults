genesOfInterest=["PPP1R2P9","CXorf61","FAM47B","MAGEB2","DDX53"]

mapGenes=mapidtoname(genesOfInterest,database=geneDatabase,idcolumn='GeneSymbol')
mapGenes.loc[:,["ProbeID","GeneSymbol"]]
geneDataPlot=geneExpressionData.loc[geneExpressionData['geneExpression.Probe_ID'].isin(mapGenes['ProbeID']),:]

""" TOTAL COHORT"""
geneDataPlot.index=geneDataPlot['geneExpression.Probe_ID']
geneData_cluster1=geneDataPlot.loc[:,clusterWithSamplesEGM['cluster1'].team].T
geneData_cluster1=geneData_cluster1.assign(Cluster="Cluster 1")
geneData_cluster2=geneDataPlot.loc[:,clusterWithSamplesEGM['cluster2'].team].T
geneData_cluster2=geneData_cluster2.assign(Cluster="Cluster 2")
geneData_cluster3=geneDataPlot.loc[:,clusterWithSamplesEGM['cluster3'].team].T
geneData_cluster3=geneData_cluster3.assign(Cluster="Cluster 3")
geneExpression_AR_combined=pd.concat([geneData_cluster1,geneData_cluster2,geneData_cluster3])
geneExpression_AR_melt=pd.melt(geneExpression_AR_combined,id_vars="Cluster")
test=geneExpression_AR_combined
del test['Cluster']
plt.figure()
sns.heatmap(test,cmap='RdYlGn_r')
plt.show()
plt.figure()
plt.title('Expression')
sns.boxplot(x='geneExpression.Probe_ID',y='value',hue='Cluster',data=geneExpression_AR_melt)
plt.xticks(rotation=90)
plt.show()

""" EPIC COHORT"""
pathToAllGeneExpression="E:/Phd Thesis/Chapter 4 subsection EGM/EGM/Genes and CpG methylation EGM/Data For ADA/data downloaded from ngs-data-2/Full EGM data set transcriptomics with AHRR.txt"
allGeneExpression=pd.read_csv(pathToAllGeneExpression,",")
allGeneExpression=pd.read_csv(pathToAllGeneExpression,",",index_col=0)
allGeneExpression=allGeneExpression.clip(lower=0)
geneExpression_AR=allGeneExpression.loc[allGeneExpression.index.isin(["A_23_P34007","A_32_P105381","A_32_P379105","A_23_P254831","A_23_P432352"]),:]
geneExpression_AR_cluster1=geneExpression_AR.loc[:,itlayDataFrame[itlayDataFrame['Cluster']=="cluster 1"].Sample].T
geneExpression_AR_cluster1=geneExpression_AR_cluster1.assign(Cluster="Cluster 1")
geneExpression_AR_cluster2=geneExpression_AR.loc[:,itlayDataFrame[itlayDataFrame['Cluster']=="cluster 2"].Sample].T
geneExpression_AR_cluster2=geneExpression_AR_cluster2.assign(Cluster="Cluster 2")
geneExpression_AR_combined=pd.concat([geneExpression_AR_cluster1,geneExpression_AR_cluster2])
geneExpression_AR_melt=pd.melt(geneExpression_AR_combined,id_vars="Cluster")
heatmap1_data = pd.pivot_table(geneExpression_AR_melt, values='value', index='Cluster', columns='Probe_ID')
plt.figure()
top50Cluster1=sns.heatmap(heatmap1_data.T,cmap='RdYlGn_r', linewidths=1.5,yticklabels=False)
top50Cluster1=sns.clustermap(heatmap1_data.T,cmap='RdYlGn_r', linewidths=1.5,yticklabels=False)
plt.show()
plt.figure()
sns.boxplot(x='Probe_ID',y='value',hue='Cluster',data=geneExpression_AR_melt)
plt.xticks(rotation=90)
plt.show()


from scipy.stats import ks_2samp
import numpy as np
ks_2samp(x, y)