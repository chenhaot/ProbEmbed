# -*- coding: utf-8 -*-
"""
Created on Sun Jul 06 19:30:32 2014

@author: Adith
"""

import nltk.corpus
import nltk.tokenize.punkt
import nltk.stem.wordnet
import nltk.tag.simplify
import collections
import string
import cPickle
import sys

stopWords = set()
#stopWFile = open('nltk_data/corpora/stopwords/english','r')
#for line in stopWFile:
#    token = line.strip().lower()
#    stopWords.add(token)
#stopWFile.close()

lmtzr = nltk.stem.wordnet.WordNetLemmatizer()
text = nltk.corpus.brown.tagged_sents()

wordNetTagMap = {'ADJ':'a', 'ADV':'r', 'FW':'n',
                 'MOD':'v', 'N':'n', 'NP':'n',
                 'NUM':'n', 'PRO':'n', 'V':'v',
                 'VD':'v', 'VG':'v', 'VN':'v'}
                 
def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False

#punctTransTable = string.maketrans(string.punctuation, ' '*len(string.punctuation))

def processWordTag(word, tag):
    if len(word) > 1:
        lower_word = word.lower()
        if lower_word in stopWords:
            return '_SSS_'
        new_tag = wordNetTagMap.get(nltk.tag.simplify.simplify_brown_tag(tag), 'n')
        new_word = lmtzr.lemmatize(lower_word, new_tag)
        final_word = new_word.translate(None, string.punctuation).strip()
        if isNumber(final_word):
            return '_NNN_'
        if len(final_word)>1:
            return new_word
    return None

trigramCounts = collections.Counter()

print len(text),
completed = 0      
for tagged_sent in text:
    simplified = [ processWordTag(word, tag) for word, tag in tagged_sent]
    simplified = filter(None, simplified)
    for trg in nltk.trigrams(simplified):
        trigramCounts[trg]+=1
    completed+=1
    if completed%100 == 0:
        print ".",
        sys.stdout.flush()
f = open("Brown_trigrams.pkl",'wb')
cPickle.dump(trigramCounts, f, -1)
f.close()
