class createProfilePlot:

    def __init__(self,name):
        self.plot_name=name
        self.data = None

    def set_data_set(self,dataset):
        self.data=dataset

    def specify_groups_for_barplot(self,groups):
        # groups can be an numpy array
        group_code=[]
        for x in range(groups.__len__()):
            for y in range(groups[x].__len__()):
                group_code.append(x)
        self.group_code=group_code
