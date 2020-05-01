import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def extractfeatures(Wmatrix, methodScore,score_threshold, clustersize):
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

def createdisplotofscoredfeatures(scoredmatrix=None,title=None,xlabel=None,ylabel=None):
    featureDistPlot=sns.distplot(scoredmatrix)
    featureDistPlot.axes.set_title(title, fontsize=16)
    featureDistPlot.set_xlabel(xlabel, fontsize=14)
    featureDistPlot.set_ylabel(ylabel, fontsize=14)
    featureDistPlot.tick_params(labelsize=10)
    plt.show()
    return featureDistPlot

def determineclustermembers(Hmatrix, samplenames):
    #For each column in H, maximum value will be determined and right cluster assigned
    index_maxvalue = np.argmax(Hmatrix, axis=0)
    ClusterMembers = []
    for cluster in range(np.min(index_maxvalue), np.max(index_maxvalue) + 1):
        ClusterMembers.append([i for indx, i in enumerate(samplenames) if index_maxvalue[indx] == cluster])
    return ClusterMembers, index_maxvalue