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


# this is unnecessary
'''
f = open('../data/Brown_vocab.pkl','rb')
tup = cPickle.load(f)
vocabulary = tup[0]
vocabList = tup[1]
f.close()
'''

numDim = 2
if len(sys.argv) > 1:
    numDim = int(sys.argv[1])

reg_lambda = 1.0
if len(sys.argv) > 2:
    reg_lambda = float(sys.argv[2])

dp_lambda = 1.0
if len(sys.argv) > 3:
    dp_lambda = float(sys.argv[3])


f = open('../data/Brown_prob','rb')
npzFiles = numpy.load(f)
true_probabilities = scipy.sparse.coo_matrix((npzFiles['arr_2'], (npzFiles['arr_0'], npzFiles['arr_1'])))
numWords = true_probabilities.shape[0]
probabilities = true_probabilities.toarray()
f.close()

ops = {'maxiter': 10000, 'disp': True, 'ftol': 1e-5, 'gtol': 1e-5, 'maxcor': 100}

def DPMeans(currEmbedding):
    global numWords, dp_lambda
    print "\n"
    sys.stdout.flush()

    clusterAssignments = numpy.zeros(numWords)
    # initial there is only one cluster
    clusterMeans = numpy.mean(currEmbedding, axis=0)
    clusterMeans = numpy.atleast_2d(clusterMeans)

    potential = None
    prev_potential = None
    improvement = 1
    # add an error tolerance
    while improvement > ops['ftol']:
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
                # how can this happen
                print 'error in dp means, no clusters can be empty'
                clusterMeans[i,:] = 0
        if potential == 0:
          break
        if prev_potential is not None:
          improvement = (prev_potential-potential)/prev_potential
        prev_potential = potential;
        print 'clustering', potential, numpy.shape(clusterMeans)[0]
        sys.stdout.flush()

    return (clusterMeans, clusterAssignments)

def KLDiv(currEmbeddingFlat, clusterMeans, clusterAssignments):
    global reg_lambda, dp_lambda, probabilities, numWords, numDim
    currEmbedding = numpy.reshape(currEmbeddingFlat, (numWords, numDim))

    distanceMatrix = 1 + scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(currEmbedding, 'sqeuclidean'))
    logDistanceMatrix = numpy.log(distanceMatrix)
    weightedDistanceMatrix = numpy.multiply(probabilities, logDistanceMatrix)
    # minimize KL divergence, this should be minus
    Objective = weightedDistanceMatrix.sum()

    invDistanceMatrix = numpy.reciprocal(distanceMatrix)
    Z = invDistanceMatrix.sum(axis=1)
    logZ = numpy.log(Z)

    # this should also be minus
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
            # a bug here, missed reg_lambda
            gradient[k,:] += 2.0*reg_lambda*(wordEmbed-clusterEmbed)

    gradientFlat = numpy.reshape(gradient, numWords*numDim)
    return (Objective, gradientFlat)


# add numDim so that dp_lambda roughly makes sense
X = numpy.random.randn(numWords, numDim)/numpy.sqrt(numDim)
# for i in xrange(10):        #TODO: Truly iterate to convergence

prevObjective, currObjective = None, None
improvement, iteration = 1, 0
while improvement > ops['ftol']:
    clusterMeans, clusterAssignments = DPMeans(X)
    f = open('../data/Brown_embed.%d.%d' % (iteration, numDim),'wb')
    numpy.savez(f, X, clusterMeans, clusterAssignments)
    f.close()

    start_x = numpy.reshape(X, numWords*numDim)
    Result = scipy.optimize.minimize(fun=KLDiv, x0=start_x,
        args=(clusterMeans, clusterAssignments), method='L-BFGS-B',
        jac=True, tol=ops['ftol'], options=ops)
    
    X = numpy.reshape(Result['x'], (numWords, numDim))
    # only 1 value
    currObjective = Result['fun'].sum()
    if prevObjective is not None:
      improvement = (prevObjective - currObjective) / prevObjective
    prevObjective = currObjective
    print "ITERATION", iteration, Result['success'], Result['status'], Result['message'], Result['nit'], improvement
    iteration += 1
    sys.stdout.flush()

clusterMeans, clusterAssignments = DPMeans(X)
f = open('../data/Brown_result.%d' % numDim,'wb')
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
