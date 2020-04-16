
# Read the file with the cluster info
clusterInfo = read.delim("E:/CCLE manusscript/Data/Results/Results used for paper/NMF clustering results/ClusterInfo_CCLE_8clusters.txt",sep='\t',stringsAsFactors=FALSE)

# Read the meta data to get the histology classification for each sample
#metaData=read.delim("E:/CCLE manusscript/Data/Results/Results used for paper/Microsatellite_Instability_Info.csv",sep=";",stringsAsFactors=FALSE)
metaData=read.delim("E:/GIT repos/visualizeNMFResults/Radial tree plot/Cell_lines_annotations_20181226.txt",sep='\t',stringsAsFactors=FALSE)
# set the names for the radial tree
headParent="CCLE"
headChild="Cluster"
clusters=c(1,2,3,4,5,6,7,8)

# set the parent and child nodes
centralNode=headParent
clusterNode=paste(headChild,clusters,sep=" ")
parentNode=paste(centralNode,clusterNode,sep=".")
# find the lineage of each cell
indexLineage=match(clusterInfo$team,metaData$CCLE_ID)
sampleLineage=metaData$Histology[indexLineage]
childNode=c()
for(x in 1:length(clusterInfo$team)){
	sample=clusterInfo$team[x]
	positionSample=match(clusterInfo$team[x],metaData$CCLE_ID)
	histology=metaData$Histology[positionSample]
	parent=parentNode[clusterInfo$cluster[x]]
	tmp=paste(parent,histology,sep='.')
	childNode=c(childNode,tmp)
}
# Get only the unique child nodes (otherwise you get to many seed nodes)
childNodeForRadial=unique(childNode)
childNodeForRadial=unique(childNode)
endNode=c()
for(x in 1:length(clusterInfo$team)){
	endNodeTmp=paste(childNode[x],x,sep='.')
	#endNodeTmp=paste(childNode,clusterInfo$team,sep='.')
	endNode=c(endNode,endNodeTmp)
}
# format to the input for the radial tree js script
centralNodeDf=data.frame('id'=centralNode,'weight'=0)
parentNodeDf=data.frame('id'=parentNode,'weight'=0)
childNodeDf=data.frame('id'=childNodeForRadial,'weight'=0)
finalNodeDf=data.frame('id'=endNode,'weight'=clusterInfo$cluster)
output=rbind(centralNodeDf,parentNodeDf,childNodeDf)

write.table(output,"E:/GIT repos/visualizeNMFResults/Radial tree plot/ClusterOutputRadial_noLastName.csv",sep=",",row.names=FALSE)