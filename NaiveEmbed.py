# -*- coding: utf-8 -*-
"""
Created on Mon Jul 07 00:12:22 2014

@author: Adith
"""
import cPickle
import scipy.sparse
import scipy.optimize
import numpy.random
import sys
import numpy
import scipy.spatial.distance
import itertools

f = open('Brown_vocab.pkl','rb')
tup = cPickle.load(f)
vocabulary = tup[0]
vocabList = tup[1]
f.close()

numDim = 2
if len(sys.argv) > 1:
    numDim = int(sys.argv[1])

reg_lambda = 1.0
if len(sys.argv) > 2:
    reg_lambda = float(sys.argv[2])

dp_lambda = 1.0
if len(sys.argv) > 3:
    dp_lambda = float(sys.argv[3])

pruneProb = 0.01
if len(sys.argv) > 4:
    pruneProb = float(sys.argv[4])


log_embedding = 0
if len(sys.argv) > 5:
    log_embedding = int(sys.argv[5])

numWords = len(vocabulary)

f = open('Brown_prob','rb')
npzFiles = numpy.load(f)
probabilities = scipy.sparse.coo_matrix((npzFiles['arr_2'], (npzFiles['arr_0'], npzFiles['arr_1'])), shape=(numWords, numWords))
f.close()

prevAssignments = None
prevMeans = None

def DPMeans(currEmbedding):
    global numWords, dp_lambda, prevAssignments, prevMeans

    print "\n"

    if prevAssignments is not None:
        clusterAssignments = prevAssignments
    else:
        clusterAssignments = numpy.zeros(numWords)

    if prevMeans is not None:
        clusterMeans = prevMeans
    else:
        clusterMeans = numpy.mean(currEmbedding, axis=0)
        clusterMeans = numpy.atleast_2d(clusterMeans)

    potential = -1.0
    prev_potential = -1.0
    while (potential < prev_potential) or (prev_potential < 0.0):
        prev_potential = potential;
        potential = 0.0;

        for i in xrange(numWords):
            distances = scipy.spatial.distance.cdist(numpy.atleast_2d(currEmbedding[i,:]), numpy.atleast_2d(clusterMeans), 'sqeuclidean')
            distances = distances.ravel()
            minDist = numpy.argmin(distances)
            if distances[minDist] > dp_lambda: #New cluster!
                clusterAssignments[i] = numpy.shape(clusterMeans)[0]
                newCluster = numpy.copy(currEmbedding[i,:])
                clusterMeans = numpy.vstack((clusterMeans, newCluster))
                potential += dp_lambda
            else:
                clusterAssignments[i] = minDist
                potential += distances[minDist]

        #Recompute clusters from assignments
        for i in xrange(numpy.shape(clusterMeans)[0]):
            cluster_size = numpy.count_nonzero(clusterAssignments == i)
            if cluster_size > 0:
                clusterMeans[i,:] = numpy.mean(currEmbedding[clusterAssignments == i,:], axis=0)
            else:
                clusterMeans[i,:] = 0

        print potential, numpy.shape(clusterMeans)[0]

    prevAssignments = clusterAssignments
    prevMeans = clusterMeans

    return (clusterMeans, clusterAssignments)

iterNumber = 0 
def KLDiv(currEmbeddingFlat):
    global reg_lambda, dp_lambda, probabilities, numWords, iterNumber
    currEmbedding = numpy.reshape(currEmbeddingFlat, (numWords, numDim))

    NLL = 0.0
    #gradient = numpy.zeros(numpy.shape(currEmbedding))

    processed = 0
    for i, j, v in itertools.izip(probabilities.row, probabilities.col, probabilities.data):
        processed+=1
        if processed%100000 == 0:
            print ".",
            sys.stdout.flush()

        if v < pruneProb:
            continue

        distanceMatrix = numpy.squeeze(scipy.spatial.distance.cdist(currEmbedding[i,:][None,:], currEmbedding, 'sqeuclidean')) + numpy.ones(numWords)
        invDistanceMatrix = numpy.reciprocal(distanceMatrix)
        Z = numpy.sum(invDistanceMatrix)
        nlogZ = -numpy.log(Z)

        NLL -= v*(-numpy.log(distanceMatrix[j]) + nlogZ)

        """
        for k in xrange(numWords):
            if (not (k == i)) and (not (k == j)):
                gradient[k,:] -= (-2.0*v*invDistanceMatrix[k]*invDistanceMatrix[k])*(currEmbedding[i,:] - currEmbedding[k,:])/Z

            if k == j:
                gradient[k,:] -= 2.0*v*invDistanceMatrix[j]*(currEmbedding[i,:] - currEmbedding[j,:])*(1.0 - invDistanceMatrix[j]/Z)

            if k == i:
                term1 = -invDistanceMatrix[j]*(currEmbedding[i,:] - currEmbedding[j,:])

                diffMatrix = currEmbedding[i,:] - currEmbedding
                prodMatrix1 = numpy.multiply(diffMatrix, invDistanceMatrix[:,None])
                prodMatrix2 = numpy.multiply(prodMatrix1, invDistanceMatrix[:,None])
                unscaledTerm2 = numpy.sum(prodMatrix2, axis=0)
                term2 = unscaledTerm2/Z

                gradient[k,:] -= 2.0*v*(term1 + term2)
        """
    #Regularization
    clusterMeans, assignments = DPMeans(currEmbedding)
    for k in xrange(numWords):
        wordEmbed = currEmbedding[k,:]
        clusterEmbed = clusterMeans[assignments[k],:]
        NLL += reg_lambda*scipy.spatial.distance.cdist(numpy.atleast_2d(wordEmbed), numpy.atleast_2d(clusterEmbed), 'sqeuclidean')
        NLL += reg_lambda*dp_lambda*numpy.shape(clusterMeans)[0]
        #gradient[k,:] += 2.0*(wordEmbed-clusterEmbed)

    #gradientFlat = numpy.reshape(gradient, numWords*numDim)
    #return (NLL, gradientFlat)

    iterNumber += 1
    print "\n"+str(iterNumber)

    if log_embedding > 0:
        f = open('Brown_embed.'+str(iterNumber),'wb')
        numpy.savez(f, currEmbedding, clusterMeans, assignments)
        f.close()

    return NLL

X = numpy.random.randn(numWords, numDim)
start_x = numpy.reshape(X, numWords*numDim)

print len(probabilities.data), numWords

ops = {'maxiter': 100000, 'disp': True}

Result = scipy.optimize.minimize(fun = KLDiv, x0 = start_x, method = 'Powell', jac = False, tol = 0.000001, options = ops)
print Result['success'], Result['status'], Result['message'], Result['nit']

final_x = numpy.reshape(Result['x'], (numWords, numDim))
final_clusters, final_assignments = DPMeans(final_x)
f = open('Brown_result','rb')
numpy.savez(f, final_x, final_clusters, final_assignments)
f.close()
