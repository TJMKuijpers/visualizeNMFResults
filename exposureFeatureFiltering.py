import pandas as pd
from univariatefilterinregression import univariatefiltering
from refcvvalidation import rfecvCalculation
# load the gene expression data
pathToGeneExpression="E:/Project Envirogenomarkers/Data/Data for analysis/Data for EGM analysis/Univariate filtered data for NMF/GeneExpressionData_NMF_ReducedObject.csv"
geneExpression=pd.read_csv(pathToGeneExpression,sep=",",index_col=0)
pathToExposure="E:/Project Envirogenomarkers/Data/Data for analysis/Data for EGM analysis/unfiltered Data/POPs-v5.1.txt"
exposureData=pd.read_csv(pathToExposure,sep="\t")
exposureData=exposureData.loc[exposureData['egm_id'].isin(geneExpression.columns),:]
pathMethylation="E:/Project Envirogenomarkers/Data/Data for analysis/Data for EGM analysis/MethylationDataObject_EGM.txt"
methylationData=pd.read_csv(pathMethylation,sep="\t",index_col=0)
#
exposureData=exposureData.set_index('egm_id')
exposureData=exposureData.reindex(index=geneExpression.columns)
exposureData=exposureData.reset_index()
exposurecolumns=['DDE']
#,'PCB-153', 'PCB-138', 'PCB-156','PCB-180', 'PCB-170']
compoundRelatedProbesGenes=dict()
#for x in exposurecolumns:
#    filter=univariatefiltering(X=geneExpression,y=exposureData.loc[:,x])
#    probes=geneExpression.index[filter.get_support()]
#    filteredData=geneExpression.loc[geneExpression.index.isin(probes),:]
#    RFEfilter=rfecvCalculation(X=filteredData,y=exposureData.loc[:,x])
#    filteredprobes=filteredData.index[RFEfilter.get_support()]
#    compoundRelatedProbesGenes['Compound '+x]=filteredprobes

compoundRelatedProbesCpG=dict()
for yx in exposurecolumns:
    filter=univariatefiltering(X=methylationData,y=exposureData.loc[:,yx])
    probes=geneExpression.index[filter.get_support()]
    filteredData=geneExpression.loc[methylationData.index.isin(probes),:]
    RFEfilter=rfecvCalculation(X=filteredData,y=exposureData.loc[:,yx])
    filteredprobes=filteredData.index[RFEfilter.get_support()]
    compoundRelatedProbesCpG['Compound '+yx]=filteredprobes