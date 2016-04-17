# COS424 Final Project

Our final project is an extension to the second assignment, in which we aim to impute the methylation values along an entire chromosome genome sequence for a sparsely-sampled tissue based off 33 samples from other tissues. The 33 given samples are densely-sampled.

## Overall methodology.

We try to fit a GMM model with off-diagonal bands as nonzero covariances. This is hard because default GMM implementations do not suffice for fitting this model on such a large dataset (each tissue sample has on the order of ~450K values).

## Installation

### Notice on `git clone`

Use `git clone --recursive ...` instead of just `git clone` so that the submodule containig the original compressed data is retained.

## Tasks

### `cycles` Environment Setup

1. Make sure you have a folder named just by your netid
2. clone the repo inside of that folder (https://github.com/PrincetonUniversity/final-proj-cos-424)

Who still has to do this?

* Jerry
* Sean

### `cycles` Environment Usage

1. When you cd inside to the cyc424 in cycles, run source py/bin/activate to get a python3 env with sklearn.
2. Make sure to use tmux to keep your computations alive
3. Make changes only in the cloned repo in your personal folder so you don't interfere with others' stuff, then push it when it's ready. Don't push broken code.

## Project tasks

1. Type up project goal
2. Do Step (0) from the previous email (that involves retrieving and parsing the data in the same way as assignment 2).
3. Do Steps (1) and (2) from the previous email.  This involves setting up an ipython notebook in the repo, and retrieving the original 33-chromosome data.
4. Find related works. Read them, provide quick summaries of the relevant points useful to our project.
5. Do Step (3) from the previous email.
6. Write the intro + background section of our paper (whole paper is ~8pages).
7. Write first part of methods section (the theory) of the paper.
8. Make a section on the dataset alone - describe its sparsity, number of chromosomes we're learning the clusters from, etc.


