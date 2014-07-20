ProbEmbed
=========

Implementing a probabilistic model for learning embeddings, using Dirichlet process priors and heavy-tailed distributions.

## Requirements:
*(compatibility checked with listed version numbers)*
* Python 2.7.3
* NLTK 2.0
* Numpy 1.8.1
* Scipy 0.14.0
* Matplotlib 1.1.1

## Scripts:

### CountTrigrams.py
This script loads the Brown corpus available in NLTK, and creates a Counter containing trigram counts in `Brown_trigrams.pkl` in the working directory.
Words are lemmatized to their WordNet stem, and numbers are normalized to `_NNN_`.
*[Not recommended:]* If an appropriate stopwords list is provided, their occurrences will be normalized to `_SSS_`.

### CreateReplacementDistribution.py
**Usage:**

`python CreateReplacementDistribution.py [pruneFreq, default:3]`

This script loads `Brown_trigrams.pkl` and estimates replacement probabilities of words using the distributional hypothesis.
It expects the NLTK English stopword list, update **line:19** to point to a non-standard location.
Tokens with unigram frequency fewer than `pruneFreq` are removed from the vocabulary, as are stop words.
The vocabulary is stored in `Brown_vocab.prob`, and the empirical replacement distribution in `Brown_prob.pkl`.

### NaiveEmbed.py
**Usage:**

`python NaiveEmbed.py [numDim, default:2] [lambda, default:1.0] [dp_lambda, default:1.0] [pruneProb, default:0.01] [log_embedding, default:0]`

* `numDim` is the dimensionality of the learned embedded space.
* `lambda` is the regularization parameter. Always set `lambda >= 0`, `lambda >> 0` => embedding underfits the empirical replacement distribution, `lambda ~ 0` promotes over-fitting.
* `dp_lambda` is a parameter used in the Dirichlet process prior. Always set `dp_lambda >= 0`, `dp_lambda >> 0` => all words embedded into same cluster, `dp_lambda ~ 0` => lots of less-separated clusters.
*Note:* `dp_lambda = 0` => `lambda = 0` mathematically.
* `pruneProb` is a computational hack to speed up training. Probabilities smaller than `pruneProb` are ignored when computing KL-divergence.
* `log_embedding` indicates whether the embedding after each step of gradient descent should be logged. If `log_embedding > 0`, intermediate embeddings will be stored in `Brown_embed.[iterationNumber]`

**TODO:** Clusters are cached between iterations to warm-start the DP-means algorithm, but is incorrect since DP-means monotonically splits existing clusters.

### plotEmbeddings.py
**Usage:**

`python plotEmbeddings.py [fileName, default:'Brown_embed.1']`

This script loads `Brown_vocab.pkl` and the embedding file, and plots the first two dimensions of the learned embeddings.
It also lists up to 20 words that belong to each cluster.

**TODO:** Print words closest to cluster center rather than by vocabularyID.

### ToDo: NoiseConstrastiveEmbed.py
