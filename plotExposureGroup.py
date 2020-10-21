import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
def plotExposureGroup(inputdata=None,columnToTake=None,clusterGroups=None):
    totaldata=pd.DataFrame()
    indexCluster=0
    for x in inputdata:
        count_data=pd.DataFrame({'Group':x[columnToTake].value_counts()})
        count_data=count_data.assign(cluster=clusterGroups[indexCluster])
        totaldata=pd.concat([totaldata,count_data])
        indexCluster=indexCluster+1

    # get the index of total data to be a column stating the exposure group
    totaldata['exposure group'] = totaldata.index
    # melt the data to long format for seaborn
    totaldata_melt=pd.melt(totaldata, id_vars=['exposure group', 'cluster'])
    print(totaldata_melt)
    plt.figure()
    categoryplot=sns.catplot(x='cluster',y='value',hue='exposure group',data=totaldata_melt,kind='bar')
    plt.show()
    return categoryplot