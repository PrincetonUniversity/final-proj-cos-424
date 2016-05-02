# COS424 Final Project

Our final project is an extension to the second assignment, in which we aim to impute the methylation values along an entire chromosome genome sequence for a sparsely-sampled tissue based off 33 samples from other tissues. The 33 given samples are densely-sampled.

## Overall methodology.

We try to fit a GMM model with off-diagonal bands as nonzero covariances. This is hard because default GMM implementations do not suffice for fitting this model on such a large dataset (each tissue sample has on the order of ~450K values).

## Directory description

A first-level `data/` directory may be formed locally with cached computation data or unzipped raw data. This is intentionally git-ignored so that the repository does not store cacheable values.

The `methylation-imputation/` directory is the Assigment 2 directory with starter code and zipped raw data. It is a submodule.

The `utils/` directory contains various factored helper functions.

The first-level of the git directory contains various `jupyter` notebooks that run the models in our paper in addition to the directories named above.

## Installation

### Dependencies

All the code uses Python 3; we assume typical machine learning libraries such as `pandas`, `numpy`, `scikit-learn`, etc. Please install them on-demand.

### Notice on `git clone`

Use `git clone --recursive ...` instead of just `git clone` so that the submodule containig the original compressed data is retained.

## Tasks

### `cycles` Environment Setup

1. Make sure you have a folder named just by your netid
2. clone the repo inside of that folder (https://github.com/PrincetonUniversity/final-proj-cos-424)

### `cycles` Environment Usage

1. When you cd inside to the cyc424 in cycles, run source py/bin/activate to get a python3 env with sklearn.
2. Make sure to use tmux to keep your computations alive
3. Make changes only in the cloned repo in your personal folder so you don't interfere with others' stuff, then push it when it's ready. Don't push broken code.
4. To use a `jupyter` notebook server hosted on cycles (but have graphical browser access on your local machine), see below section.

#### Remote `jupyter` notebook setup

Let `PORT` be some unused socket number (choose something between 10000 and 15000).

0. `ssh` into `cycles` with the `ssh` option `-L8888:localhost:$PORT`
1. Set up your `cycles` environment (including activating the virtual Python environment).
2. In the directory containing the notebook you wish to edit, run `jupyter notebook --no-browser --port=$PORT`
3. Open your browser locally, and navigate to `localhost:8888`

## Project tasks

1. [SEAN] Retrieve and parse data + get more data from Barbra (WGBS methylation data) - describe its sparsity, number of chromosomes we're learning the clusters from, etc.
2. [SEAN] Re-run exploratory analysis on the new data, write up section on the dataset.
3. [VLAD] Write code for diagonal GMM MLE predictor and argmax cluster predictor (also expectation predictor?)
4. [JERRY] Write TensorFlow EM GMM code for tridiagonal matrix likelihood (inner product).
5. [SID] Write methods section for GMM EM code and posterior marginalization (theory, no proofs needed).