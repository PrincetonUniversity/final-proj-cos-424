# COS424 Final Project

Our final project is an extension to the second assignment, in which we aim to impute the methylation values along an entire chromosome genome sequence for a sparsely-sampled tissue based of 34 samples from other tissues. The 34 given samples are densely-sampled.

4 of the 34 samples were pathological, and thus not used in the model.

## Overall methodology.

We try to fit a GMM model with off-diagonal bands as nonzero covariances. This is hard because default GMM implementations do not suffice for fitting this model on such a large dataset (each tissue sample has on the order of ~380K values).

We moved away from the off-diagonal model, relying only on the diagonal bands instead. See the paper for justification of the decision.

## Directory description

A first-level `data/` directory may be formed locally with cached computation data or unzipped raw data. This is intentionally git-ignored so that the repository does not store cacheable values.

The `methylation-imputation/` directory is the Assigment 2 directory with starter code and zipped raw data. It is a submodule.

The `utils/` directory contains various factored helper functions.

The `paper.pdf` file contains the submitted final paper.

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
