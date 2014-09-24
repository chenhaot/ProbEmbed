from matplotlib import pyplot as plt
import sys
import numpy
import cPickle
import scipy.sparse
import scipy.spatial.distance

f = open('Brown_vocab.pkl','rb')
(vocabulary, vocabList, contextWordDict, wordContextDict) = cPickle.load(f)
f.close()

numWords = len(vocabulary)
print numWords

f = open('Brown_prob','rb')
npzFiles = numpy.load(f)
true_probabilities = scipy.sparse.coo_matrix((npzFiles['arr_2'], (npzFiles['arr_0'], npzFiles['arr_1'])), shape=(numWords, numWords))
probabilities = true_probabilities.toarray()
f.close()

def ProbableReplaceablePairs(prob):
    global vocabList
    numpy.fill_diagonal(prob, 0)

    partitionIndices = numpy.argpartition(prob, -50, axis=None)
    partitionIndices = partitionIndices[-50:]
    unraveledIndices =  numpy.unravel_index(partitionIndices, numpy.shape(prob))
    for i in xrange(50):
        rowID = unraveledIndices[0][i]
        colID = unraveledIndices[1][i]
        print vocabList[rowID], vocabList[colID], prob[rowID, colID]

def deriveProbabilities(embed):
    distanceMatrix = 1 + scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(embed, 'sqeuclidean'))
    invDistanceMatrix = numpy.reciprocal(distanceMatrix)
    Z = invDistanceMatrix.sum(axis=1)
    currentDistribution = invDistanceMatrix / Z[:,numpy.newaxis]
    return currentDistribution

def ProbFromFile(fileName):
    f = open(fileName,'rb')
    npzFiles = numpy.load(f)
    embedding = npzFiles['arr_0']
    clusterMeans = npzFiles['arr_1']
    clusterAssignments = npzFiles['arr_2']
    f.close()

    return deriveProbabilities(embedding)

print "EMPIRICAL DISTRIBUTION"
ProbableReplaceablePairs(probabilities)
"""
print "RANDOM EMBEDDING"
randomEmbedProb = ProbFromFile('Brown_embed.1')
ProbableReplaceablePairs(randomEmbedProb)
"""
print "DPMEANS EMBEDDING"
dpEmbedProb = ProbFromFile('Brown_result')
ProbableReplaceablePairs(dpEmbedProb)
