import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
def createprofileplot(inputData=None,title=None,xlabel= None, ylabel= None):
    boxplot=sns.boxplot(data=inputData,palette="colorblind",width=0.5)
    boxplot.axes.set_title(title, fontsize=16)
    boxplot.set_xlabel(xlabel,fontsize=14)
    boxplot.set_xticklabels(boxplot.get_xticklabels(), rotation=45)
    boxplot.set_ylabel(ylabel,fontsize=14)
    boxplot.tick_params(labelsize=10)
    plt.show()

def creategroupprofileplot(inputData=None,xVar=None,yVar=None,hueVar=None,title=None,xlabel=None,ylabel=None):
    boxplotGroup=sns.boxplot(x=xVar, y=yVar, hue=hueVar, data=inputData)
    boxplotGroup.axes.set_title(title, fontsize=16)
    boxplotGroup.set_xlabel(xlabel, fontsize=14)
    boxplotGroup.set_ylabel(ylabel, fontsize=14)
    boxplotGroup.tick_params(labelsize=10)
    plt.show()

def createboxplotwithpvaluesbetweengroups(inputData=None,xVar=None,yVar=None,hueVar=None,title=None,xlabel=None,ylabel=None,significanceColumns=None,columnForMax=None):
    boxplotGroup=sns.boxplot(x=xVar, y=yVar, hue=hueVar, data=inputData)
    boxplotGroup.axes.set_title(title, fontsize=16)
    boxplotGroup.set_xlabel(xlabel, fontsize=14)
    boxplotGroup.set_ylabel(ylabel, fontsize=14)
    boxplotGroup.tick_params(labelsize=10)
    # add the significance bars between groups
    y,col=inputData[columnForMax].max(),'k'
    n=0
    for levels in significanceColumns:
        x1,x2 = levels[0],levels[1]
        h=(y/100)+(y/100)*n+(y/1000)
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        plt.text((x1 + x2) * .5, y + h, "*", ha='center', va='bottom', color=col)
        n=n+1
    plt.show()


def createcategoricalplot(inputData=None,xVar=None,yVar=None,hueVar=None,title=None):
    return None
