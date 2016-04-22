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
3. Open your browser locally, and navigate to `localhost:$PORT`

## Project tasks

1. Type up project goal
2. Retrieve and parse data + get more data from Barbra (WGBS methylation data) (Sean todo).
3. Loot at the n x n correlation matrix to get the number of clusters k.  
   Look at the correlation matrices inside the clusters. 
   This involves setting up an ipython notebook in the repo, and retrieving the original 33-chromosome data.
4. Find related works. Read them, provide quick summaries of the relevant points useful to our project.
5. Write TensorFlow EM GMM code for k clusters.
6. Write the intro + background section of our paper (whole paper is ~8pages).
7. Make a section on the dataset alone - describe its sparsity, number of chromosomes we're learning the clusters from, etc.

By 4/23:

Vlad: [3] (finish by EOD 4/17), half of [4], [7]
Sid: half of [4]
Jerry: get started on [5]
Sid and Jerry: [1] and [6]
Sean: [2]
