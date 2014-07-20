# -*- coding: utf-8 -*-
"""
Created on Mon Jul 07 00:12:22 2014

@author: Adith
"""
import cPickle
import scipy.sparse
import sys
import numpy
import scipy
import collections

f = open('Brown_trigrams.pkl','rb')
trigrams = cPickle.load(f)
f.close()

stopWords = set()
stopWFile = open('/home/adith/nltk_data/corpora/stopwords/english','r')
for line in stopWFile:
    token = line.strip().lower()
    stopWords.add(token)
stopWFile.close()

unigramFreq = collections.Counter()
for elem, cnt in trigrams.iteritems():
    prevWord = elem[0]
    currWord = elem[1]
    nextWord = elem[2]

    if currWord in stopWords:
        continue
    unigramFreq[currWord] += cnt

pruneFreq = 3
if len(sys.argv) > 1:
    pruneFreq = int(sys.argv[1])

#Prune vocabulary, only >3
toRemove = set()
for word, cnt in unigramFreq.iteritems():
    if cnt <= pruneFreq:
        toRemove.add(word)

print len(unigramFreq), len(toRemove)
vocabulary = {}
vocabList = []
wordContextDict = {}
contextWordDict = {}
for elem, cnt in trigrams.items():
    #if cnt <= 1:
    #    continue
    prevWord = elem[0]
    currWord = elem[1]
    nextWord = elem[2]

    if currWord in stopWords:
        continue
    if currWord in toRemove:
        continue

    #if prevWord not in vocabulary:
    #    vocabulary[prevWord] = len(vocabulary)
    #    vocabList.append(prevWord)
    if currWord not in vocabulary:
        vocabulary[currWord] = len(vocabulary)
        vocabList.append(currWord)
    #if nextWord not in vocabulary:
    #    vocabulary[nextWord] = len(vocabulary)
    #    vocabList.append(nextWord)

    if currWord not in wordContextDict:
        wordContextDict[currWord] = {'_SUM_':0}
    context = (prevWord, nextWord)
    if context not in contextWordDict:
        contextWordDict[context] = {'_SUM_':0}
        
    currDict = wordContextDict[currWord]
    currDict[context] = cnt
    currDict['_SUM_'] += cnt
    wordContextDict[currWord] = currDict
    
    currDict = contextWordDict[context]
    currDict[currWord] = cnt
    currDict['_SUM_'] += cnt
    contextWordDict[context] = currDict

vocab_size = len(vocabulary)
print "Vocabulary", vocab_size

possiblePairs = set()
for context, val in contextWordDict.iteritems():
    for k1 in val:
        if k1 == '_SUM_':
            continue
        i = vocabulary[k1]
        for k2 in val:
            if k2 == '_SUM_':
                continue
            j = vocabulary[k2]
            if (i,j) not in possiblePairs:
                possiblePairs.add((i,j))

print len(possiblePairs)
probabilities = scipy.sparse.dok_matrix((vocab_size, vocab_size))
processed = 0
for i,j in possiblePairs:
    currProb = 0.0
    source_word = vocabList[i]
    source_dict = wordContextDict[source_word]
    target_word = vocabList[j]
    for key, val in source_dict.iteritems():
        if key == '_SUM_':
            continue
        if target_word not in contextWordDict[key]:
            continue
        tempDict = contextWordDict[key]
        currProb += (val*1.0/source_dict['_SUM_'])*(tempDict[target_word]*1.0/tempDict['_SUM_'])
    if currProb > 0.0:
        probabilities[i, j] = currProb
    processed += 1
    if processed%1000 == 0:
        print ".",
        sys.stdout.flush()

f = open('Brown_vocab.pkl','wb')
cPickle.dump((vocabulary, vocabList, contextWordDict, wordContextDict), f, -1)
f.close()

probabilities = probabilities.tocoo()
f = open('Brown_prob','wb')
numpy.savez(f, probabilities.row, probabilities.col, probabilities.data)
f.close()
