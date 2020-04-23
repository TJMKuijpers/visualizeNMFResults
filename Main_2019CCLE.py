import pandas as pd
import numpy as np
import pickle
import seaborn as sns
from getMetaDataInfo import getMetaDataInfo

def extractFeatures(Wmatrix, methodScore,score_threshold, clustersize):
	# Script to score the basis matrix
	print('Get the basis matrix')
	scores = np.zeros(Wmatrix.shape[0])
	print('score features')
	if (methodScore == 'Kim'):
		# Part to score the matrix
		for x in range(Wmatrix.shape[0]):
			probability = Wmatrix[x, :] / Wmatrix[x, :].sum()
			scores[x] = np.dot(probability, np.log2(probability + 0.0001).T)
		scores = 1. + 1. / np.log2(Wmatrix.shape[1]) * scores
		scoredMatrix = scores
		threshold = np.median(scores) + 3 * np.median(abs(scores - np.median(scores)))
		sel = scores > threshold
		m = np.median(scores)
		sel = np.array([sel[i] and np.max(Wmatrix[i, :]) > m for i in range(Wmatrix.shape[0])])
		features_per_cluster = [[] for i in range(Wmatrix.shape[1])]
		for i in range(Wmatrix.shape[0]):
			if sel[i] == True:
				pos = np.argmax(Wmatrix[i, :])
				features_per_cluster[pos].append(i)

	if (methodScore == 'Max_Kim_Score'):
		scoredMatrix = np.argmax(Wmatrix, axis=1)
		# Part to score the matrix
		for x in range(Wmatrix.shape[0]):
			probability = Wmatrix[x, :] / Wmatrix[x, :].sum()
			scores[x] = np.dot(probability, np.log2(probability + 0.0001).T)
		scores = 1. + 1. / np.log2(Wmatrix.shape[1]) * scores
		scoredMatrix = scores
		sel = scores >=score_threshold
		m = np.median(scores)
		sel = np.array([sel[i] and np.max(Wmatrix[i, :]) > m for i in range(Wmatrix.shape[0])])
		features_per_cluster = [[] for i in range(Wmatrix.shape[1])]
		for i in range(Wmatrix.shape[0]):
			if sel[i] == True:
				pos = np.argmax(Wmatrix[i, :])
				features_per_cluster[pos].append(i)
	return scoredMatrix, sel, features_per_cluster

def determineClusterMembers(Hmatrix, samplenames):
	"For each column in H, maximum value will be determined and right cluster assigned"
	index_maxvalue = np.argmax(Hmatrix, axis=0)
	ClusterMembers = []
	for cluster in range(np.min(index_maxvalue), np.max(index_maxvalue) + 1):
		ClusterMembers.append([i for indx, i in enumerate(samplenames) if index_maxvalue[indx] == cluster])
	return ClusterMembers, index_maxvalue

# First read the file that contains the list of features for each cluster
# Note that for the paper only cluster 6 and cluster 7 are analyzed
path_to_genes_cluster6 = "E:/CCLE manusscript/Data\Results/Results used for paper/NMF clustering results/Genes_cluster6_k8_nrun40.csv"
cluster6_feature_genes = pd.read_csv(path_to_genes_cluster6,sep=',',index_col=False)
path_to_genes_cluster7 = "E:/CCLE manusscript/Data/Results/Results used for paper/NMF clustering results/Genes_cluster7_k8_nrun40.csv"
cluster7_feature_genes = pd.read_csv(path_to_genes_cluster7,sep=",",index_col=False)

# Read the gene expression data
path_to_gene_expression_file = "E:/CCLE manusscript/Data/Cluster results/GeneExpressionSet_CCLE.csv"
gene_expression_data=pd.read_csv(path_to_gene_expression_file,sep=",",index_col=0)
# Load the W matrix to score the features with only features equal to one
pickle_in = open("E:/CCLE manusscript/Data/Results/Results used for paper/NMF clustering results/Wmatric_40runs_pickle","rb")
Wmatrix=pickle.load(pickle_in)
pickle_in1=open("E:/CCLE manusscript/Data/Results/Results used for paper/NMF clustering results/Hmatrix_40runs.pickle","rb")
Hmatrix=pickle.load(pickle_in1)
path_to_gene_expression_file = "E:/CCLE manusscript/Data/Cluster results/GeneExpressionSet_CCLE.csv"
gene_expression_data=pd.read_csv(path_to_gene_expression_file,sep=",",index_col=0)
# Load the methylation data
path_to_methylation_file="E:/CCLE manusscript/Data/Cluster results/MethylationDataSet_CCLE.csv"
methylationdata=pd.read_csv(path_to_methylation_file,sep=',',index_col=0)
# Load the meta data
file_microsatellite="E:/CCLE manusscript/Data/Results/Results used for paper/Microsatellite_Instability_Info.csv"
microsatellite=pd.read_csv(file_microsatellite,sep=";")
meta_data_lineage="E:/CCLE manusscript/Data/RAW data/Cell_lines_annotations_20181226.txt"
meta_data_lineage=pd.read_csv(meta_data_lineage,sep='\t')

# The Wmatrix and H matrix with the lowest error are determined from the NMF results earlier
cpg_wmatrix=Wmatrix[35][1]
genes_wmatrix=Wmatrix[35][0]
scored_matrix_genes,sel_genes,features_per_cluster_genes=extractFeatures(genes_wmatrix,'Kim',score_threshold=0.95,clustersize=8)
scored_matrix_cpg,sel_cpg,features_per_cluster_cpg=extractFeatures(cpg_wmatrix,'Kim',score_threshold=0.8,clustersize=8)
scored_matrix_maximum_genes,sel_genes,features_maximum_per_cluster_genes=extractFeatures(genes_wmatrix,'Max_Kim_Score',score_threshold=0.95,clustersize=8)
scored_matrix_maximum_cpg,sel_cpg,features_maximum_per_cluster_cpg=extractFeatures(cpg_wmatrix,'Max_Kim_Score',score_threshold=0.75,clustersize=8)
hmatrix_minimal_error=Hmatrix[35]
cluster_members,index_maxvalue=determineClusterMembers(hmatrix_minimal_error,samplenames=gene_expression_data.columns.values)




analyze_haematopoietic_lymphoid_tissue=getMetaDataInfo(cluster_number='cluster7')
analyze_haematopoietic_lymphoid_tissue.samples_in_cluster(cluster_members[6])
analyze_haematopoietic_lymphoid_tissue.set_meta_data(meta_data=microsatellite)
analyze_haematopoietic_lymphoid_tissue.extractMeta(column='CCLE_ID')


analyze_cluster_3=getMetaDataInfo(cluster_number='cluster3')
analyze_cluster_3.samples_in_cluster(cluster_members[2])
analyze_cluster_3.set_meta_data(meta_data=microsatellite)
analyze_cluster_3.extractMeta(column='CCLE_ID')
analyze_cluster_4=getMetaDataInfo(cluster_number='cluster4')
analyze_cluster_4.samples_in_cluster(cluster_members[3])
analyze_cluster_4.set_meta_data(meta_data=microsatellite)
analyze_cluster_4.extractMeta(column='CCLE_ID')
analyze_cluster_5=getMetaDataInfo(cluster_number='cluster5')
analyze_cluster_5.samples_in_cluster(cluster_members[4])
analyze_cluster_5.set_meta_data(meta_data=microsatellite)
analyze_cluster_5.extractMeta(column='CCLE_ID')

testdf=methylationdata.iloc[features_maximum_per_cluster_cpg[6],:]
test2=gene_expression_data.iloc[features_maximum_per_cluster_genes[6],:]

methlation_cluster5=methylationdata.iloc[features_maximum_per_cluster_cpg[4],:]
gene_cluster5=gene_expression_data.iloc[features_maximum_per_cluster_genes[4],:]
