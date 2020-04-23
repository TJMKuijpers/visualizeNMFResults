import seaborn as sns
import matplotlib.pyplot as plt

def createprofileplot(inputData=None,title=None,xlabel= None, ylabel= None):
    boxplot=sns.boxplot(data=inputData,palette="colorblind",width=0.5)
    boxplot.axes.set_title(title, fontsize=16)
    boxplot.set_xlabel(xlabel,fontsize=14)
    boxplot.set_xticklabels(boxplot.get_xticklabels(), rotation=45)
    boxplot.set_ylabel(ylabel,fontsize=14)
    boxplot.tick_params(labelsize=10)
    plt.show()

def creategroupprofileplot(inputData=None,title=None,xlabel=None,ylabel=None):
    boxplotGroup=sns.boxplot(data=inputData,palette="coorblind",width=0.5)
    boxplotGroup.axes.set_title(title, fontsize=16)
    boxplotGroup.set_xlabel(xlabel, fontsize=14)
    boxplotGroup.set_ylabel(ylabel, fontsize=14)
    boxplotGroup.tick_params(labelsize=10)
    plt.show()
