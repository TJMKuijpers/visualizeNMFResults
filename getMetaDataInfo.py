class getMetaDataInfo:

    def __init__(self,cluster_number):
        self.cluster_number=cluster_number
        self.meta_data = None
        self.samples = None
        self.data_for_cluster = None

    def set_meta_data(self,meta_data):
        self.meta_data=meta_data

    def samples_in_cluster(self,cluster_content):
        self.samples=cluster_content

    def extractMeta(self,column):
        subdata = self.meta_data.loc[self.meta_data[column].isin(self.samples),:]
        self.data_for_cluster=subdata