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

numWords = len(vocabulary)

f = open('Brown_prob','rb')
npzFiles = numpy.load(f)
true_probabilities = scipy.sparse.coo_matrix((npzFiles['arr_2'], (npzFiles['arr_0'], npzFiles['arr_1'])), shape=(numWords, numWords))
probabilities = true_probabilities.toarray()
f.close()

def DPMeans(currEmbedding):
    global numWords, dp_lambda

    print "\n"
    sys.stdout.flush()

    clusterAssignments = numpy.zeros(numWords)

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
        sys.stdout.flush()

    return (clusterMeans, clusterAssignments)

iterNumber = 0 
def KLDiv(currEmbeddingFlat, clusterMeans, clusterAssignments):
    global reg_lambda, dp_lambda, probabilities, numWords, numDim, iterNumber
    currEmbedding = numpy.reshape(currEmbeddingFlat, (numWords, numDim))

    distanceMatrix = 1 + scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(currEmbedding, 'sqeuclidean'))
    logDistanceMatrix = numpy.log(distanceMatrix)
    weightedDistanceMatrix = numpy.multiply(probabilities, logDistanceMatrix)
    Objective = weightedDistanceMatrix.sum()

    invDistanceMatrix = numpy.reciprocal(distanceMatrix)
    Z = invDistanceMatrix.sum(axis=1)
    logZ = numpy.log(Z)

    Objective += logZ.sum()

    currentDistribution = invDistanceMatrix / Z[:,numpy.newaxis]
    gradient = numpy.zeros(numpy.shape(currEmbedding))
    for k in xrange(numWords):
        deltaEmbedding = 2*(currEmbedding[k,:] - currEmbedding)
        firstTermProb = probabilities[:,k] - currentDistribution[:,k]
        secondTermProb = probabilities[k,:] - currentDistribution[k,:]
        deltaProb = firstTermProb + secondTermProb
        weightedProb = numpy.multiply(deltaProb, invDistanceMatrix[k,:])
        perInstanceGradient = numpy.multiply(deltaEmbedding, weightedProb[:,numpy.newaxis])
        gradient[k,:] = perInstanceGradient.sum(axis=0)

    #Regularization
    if reg_lambda > 0:
        for k in xrange(numWords):
            wordEmbed = currEmbedding[k,:]
            clusterEmbed = clusterMeans[clusterAssignments[k],:]
            Objective += reg_lambda*scipy.spatial.distance.cdist(numpy.atleast_2d(wordEmbed), numpy.atleast_2d(clusterEmbed), 'sqeuclidean')
            Objective += reg_lambda*dp_lambda*numpy.shape(clusterMeans)[0]
            gradient[k,:] += 2.0*(wordEmbed-clusterEmbed)

    gradientFlat = numpy.reshape(gradient, numWords*numDim)
    return (Objective, gradientFlat)

ops = {'maxiter': 100000, 'disp': False, 'ftol': 1e-12, 'gtol': 1e-12, 'maxcor': 100}

X = numpy.random.randn(numWords, numDim)

for i in xrange(10):        #TODO: Truly iterate to convergence
    clusterMeans, clusterAssignments = DPMeans(X)
    f = open('Brown_embed.'+str(i),'wb')
    numpy.savez(f, X, clusterMeans, clusterAssignments)
    f.close()

    start_x = numpy.reshape(X, numWords*numDim)
    Result = scipy.optimize.minimize(fun = KLDiv, x0 = start_x, args=(clusterMeans, clusterAssignments), method = 'L-BFGS-B', jac = True, tol = 1e-12, options = ops)
    print "ITERATION", i, Result['success'], Result['status'], Result['message'], Result['nit']
    sys.stdout.flush()

    X = numpy.reshape(Result['x'], (numWords, numDim))

clusterMeans, clusterAssignments = DPMeans(X)
f = open('Brown_result','wb')
numpy.savez(f, X, clusterMeans, clusterAssignments)
f.close()


"""
def ObjFun(w):
    global startClusters, startAssignments
    retVal = KLDiv(w, startClusters, startAssignments)
    return retVal[0]

def GradientFun(w):
    global startClusters, startAssignments
    retVal = KLDiv(w, startClusters, startAssignments)
    return retVal[1]

for i in xrange(5):
    startParam = numpy.random.randn(numWords*numDim)
    print scipy.optimize.check_grad(ObjFun, GradientFun, startParam)
"""
