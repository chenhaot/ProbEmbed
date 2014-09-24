from matplotlib import pyplot as plt
import sys
import numpy
import cPickle

fileName = 'Brown_embed.0'
if len(sys.argv) > 1:
    fileName = sys.argv[1]

f = open(fileName,'rb')
npzFiles = numpy.load(f)
embedding = npzFiles['arr_0']
clusterMeans = npzFiles['arr_1']
clusterAssignments = npzFiles['arr_2']
f.close()

f = open('Brown_vocab.pkl','rb')
(vocabulary, vocabList, contextWordDict, wordContextDict) = cPickle.load(f)
f.close()

colors =['b', 'g', 'r', 'c', 'm', 'y', 'k', 'brown', 'blueviolet', 'cadetblue', 'aquamarine', 'chartreuse', 'coral', 'crimson', 'cornflowerblue', 'darkcyan', 'goldenrod', 'hotpink']
#colors =['b']

numClusters = numpy.shape(clusterMeans)[0]
print "NumClusters ", numClusters
for i in xrange(numClusters):
    colorID = i%len(colors)
    plt.scatter(clusterMeans[i,0], clusterMeans[i,1], marker = 'p', edgecolors = 'none', s=10.0, color = colors[colorID])
    currClusterPoints = (clusterAssignments == i)
    currentPoints = embedding[currClusterPoints, :]
    print "Cluster ", i, numpy.shape(currentPoints)[0]
    words = numpy.flatnonzero(currClusterPoints)
    for j in xrange(min(20, len(words))):
        print vocabList[words[j]],
    print "\n"
    plt.scatter(currentPoints[:,0], currentPoints[:,1], marker = 'o', edgecolors = 'none', s=1.0, color = colors[colorID])

plt.show()
