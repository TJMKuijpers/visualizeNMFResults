################# Python script for statistical significant difference boxplot ######################
#  Author: T.J.M. Kuijpers data: 29-April-2020
from statsmodels.graphics.gofplots import qqplot
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import shapiro,mannwhitneyu

def calculatesignificantdifference():
    return None

def reportmeanandstandarddeviation(inputData):
    mean=np.mean(inputData)
    std=np.std(inputData)
    meanAndStandardDeviation=[mean,std]
    return meanAndStandardDeviation

def testnormaldistribution(inputData=None):
    # Shapiro-wilk test
    stat, p = shapiro(inputData)
    shapiroWilkResults='Statistics= '+str(stat)+"p= "+str(p)
    return shapiroWilkResults

def nonparamtericmannwhitneyutest(dataGroup1=None,dataGroup2=None):
    statistic,pvalue=mannwhitneyu(dataGroup1,dataGroup2)
    mannwhitneyuResults="Statistics= "+str(statistic)+" pvalue= "+str(pvalue)
    return mannwhitneyuResults

def qqplotfordistribution(inputData):
    qqplotFigure=qqplot(inputData, line='s')
    plt.show()
    return qqplotFigure