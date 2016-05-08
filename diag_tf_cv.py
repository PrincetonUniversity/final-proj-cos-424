#get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import tensorflow as tf
import numpy as np
import pandas as pd
from itertools import *
import sklearn
import math
import random
import sys
import multiprocessing
import scipy
from joblib import Parallel, delayed
import threading
nproc = max(1, multiprocessing.cpu_count() - 1)

if 'utils' not in sys.path:
    sys.path.append('utils')

import data_loader

# Warnings

import warnings
warnings.filterwarnings('ignore')

# Idempotent, cached data retrieval script
#print(data_loader.load_chromosome.__doc__)
train_df, test_df, train_ix, test_ix, train_tissues, tfs =     data_loader.load_chromosome_cached('1')
    
sess = tf.InteractiveSession()

# process input
if len(sys.argv) < 5: 
  raise Exception("Correct format: python diag_tf_cv.py <K> <'mle'/'argmax'> <MINCOV> <OUTPUT_FILE>")

K = int(sys.argv[1])
IS_MLE = False
if sys.argv[2] == "mle": 
  IS_MLE = True
MIN_COVAR = float(sys.argv[3]) 
OUTPUT_FILE = sys.argv[4]


# In[3]:

# Perhaps there are obvious sequence trends?
def local_impute(data):
    #http://stackoverflow.com/questions/9537543/replace-nans-in-numpy-array-with-closest-non-nan-value
    mask = np.isnan(data)
    data[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])
    return data

# Do mean imputation on our training data.
def mean_impute(data):
    mask = np.isnan(data)
    data[mask] = float(data.mean()) # just = m messes with serialization
    return data
train_df_imp = train_df
train_df_int = train_df
for i in train_tissues:
    train_df_imp[i] = mean_impute(train_df[i].copy())
    train_df_int[i] = local_impute(train_df[i].copy())


# In[4]:

# Copied pretty much directly from sklearn

# Returns a TensorFlow scalar with the size of the i-th dimension for
# the parameter tensor x.
def tf_get_shape(x, i):
    return tf.squeeze(tf.slice(tf.shape(x), [i], [1])) 

def tf_nrows(x):
    return tf_get_shape(x, 0)

def tf_ncols(x):
    return tf_get_shape(x, 1)

# Simultaneous K-cluster likelihood computation.
# X is NxD, mus is KxD, sigmas is KxD
# Output is KxN likelihoods for each sample in each cluster.
def tf_log_normals(X, mus, sigmas):
    # p(X) = sqrt(a * b * c)
    # a = (2 pi)^(-p)
    # b = det(sigma)^(-1)
    # c = exp(-(x - mu)^T sigma^(-1) (x - mu)) [expanded for numerical stability]
    #
    # Below we make simplifications since sigma is diag
    
    D = tf_ncols(mus)
    XT = tf.transpose(X) # pxN
    invsig = tf.inv(sigmas)
    
    loga = -tf.cast(D, 'float64') * tf.log(tf.constant(2 * np.pi, dtype='float64')) # scalar
    logb = tf.reduce_sum(tf.log(invsig), 1, keep_dims=True) # Kx1
    logc =          - tf.reduce_sum(invsig * tf.square(mus), 1, keep_dims=True)         + 2 * tf.matmul(invsig * mus, XT)         - tf.matmul(invsig, tf.square(XT)) # KxN
    
    return 0.5 * (loga + logb + logc)

# Stably log-sum-exps likelihood along rows.
# Reduces KxN tensor L to 1xN tensor
def tf_log_sum_exp(L):
    maxs = tf.reduce_max(L, 0, keep_dims=True) # 1xN
    return tf.log(tf.reduce_sum(tf.exp(L - maxs), 0, keep_dims=True)) + maxs

# X is NxD, mus is KxD, sigmas KxD, alphas is K
# output is KxN log likelihoods.
def tf_log_likelihood(X, mus, sigmas, alphas):
    alphas = tf.expand_dims(alphas, 1) # Kx1
    return tf_log_normals(X, mus, sigmas) + tf.log(alphas) # KxN

# X is NxD, mus is KxD, sigmas KxD, alphas is K
# output is 1xN log probability for each sample, KxN responsibilities
def estep(X, mus, sigmas, alphas):
    log_likelihoods = tf_log_likelihood(X, mus, sigmas, alphas)
    sample_log_prob = tf_log_sum_exp(log_likelihoods) # 1xN
    return sample_log_prob, tf.exp(log_likelihoods - sample_log_prob)

EPS = np.finfo(float).eps
#MIN_COVAR = 0.1
# X is NxD, resp is KxN (and normalized along axis 0)
# Returns maximize step means, covariance, and cluster priors,
# which have dimension KxD, KxD, and K, respectively
def mstep(X, resp):
    weights = tf.reduce_sum(resp, 1) # K
    invweights = tf.expand_dims(tf.inv(weights + 10 * EPS), 1) # Kx1
    alphas = EPS + weights / (tf.reduce_sum(weights) + 10 * EPS) # K
    weighted_cluster_sum = tf.matmul(resp, X) # KxD 
    mus = weighted_cluster_sum * invweights
    avg_X2 = tf.matmul(resp, tf.square(X)) * invweights
    avg_mu2 = tf.square(mus)
    avg_X_mu = mus * weighted_cluster_sum * invweights
    sigmas = avg_X2 - 2 * avg_X_mu + avg_mu2 + MIN_COVAR
    # (x - mu) (x-mu)^T for banded. 
    return mus, sigmas, alphas


def normal_likelihoods(X, mu, sigma):
    exponent = -np.dot((X - mu[np.newaxis, :]) ** 2, 1 / sigma) / 2
    return (2 * np.pi) ** (-D / 2) * np.prod(sigma) ** (-1 / 2) * np.exp(exponent)

# In[6]:

# Similar pattern to
# https://gist.github.com/narphorium/d06b7ed234287e319f18

#todo try initializing covar to emprical cv computed from kmeans labels

# Runs up to max_steps EM iterations, stopping earlier if log likelihood improves
# less than tol.
# X should be an NxD data matrix, initial_mus should be KxD
# max_steps should be an int, tol should be a float.
def fit_em(X, initial_mus, max_steps, tol, sess):
    N, D = X.shape
    K, Dmu = initial_mus.shape
    assert D == Dmu
        
    mus0 = initial_mus
    sigmas0 = np.tile(np.var(X, axis=0), (K, 1))
    alphas0 = np.ones(K) / K
    X = tf.constant(X)
    
    mus, sigmas, alphas = (tf.Variable(x, dtype='float64') for x in [mus0, sigmas0, alphas0])
    
    sess.run(tf.initialize_all_variables())
    all_ll, resp = estep(X, mus, sigmas, alphas)
    cmus, csigmas, calphas = mstep(X, resp)
    update_mus_step = tf.assign(mus, cmus)
    update_sigmas_step = tf.assign(sigmas, csigmas)
    update_alphas_step = tf.assign(alphas, calphas)     
    
    init_op = tf.initialize_all_variables()
    ll = prev_ll = -np.inf
                         
    with tf.Session() as sess2:
        sess2.run(init_op)
        for i in range(max_steps):
            ll = sess2.run(tf.reduce_mean(all_ll))
            sess2.run((update_mus_step, update_sigmas_step, update_alphas_step))
#             print('EM iteration', i, 'log likelihood', ll)
            if abs(ll - prev_ll) < tol:
                break
            prev_ll = ll
        m, s, a = sess2.run((mus, sigmas, alphas))
    
    return ll, m, s, a


# In[7]:

# Given a set of partial observations xs each of dimension O < D for a fitted GMM model with 
# K cluster priors alpha, KxD means mus, and KxD diagonal covariances sigmas,
# returns the weighted sum of normals for the remaining D - O dimensions.
#
# Returns posterior_mus, posterior_sigmas, posterior_prior,
# of dimensions:
# Kx(D-O), Kx(D-O), NxK, respectively (each mu, sigma is the same for all posteriors).
# NxK, NxKxD, NxKxD, respectively, for each x in xs, total of N.
def marginal_posterior(xs, mus, sigmas, alphas, sess):
    # https://gbhqed.wordpress.com/2010/02/21/conditional-and-marginal-distributions-of-a-multivariate-gaussian/
    # diagonal case is easy:
    # https://en.wikipedia.org/wiki/Schur_complement#Applications_to_probability_theory_and_statistics
    O = xs.shape[1]
    D = mus.shape[1]
    observed_mus, observed_sigmas = (tf.constant(a, dtype='float64')
                                     for a in (mus[:,0:O], sigmas[:, 0:O]))
    ll = tf_log_likelihood(xs, observed_mus, observed_sigmas, alphas) # KxN
    norm = tf_log_sum_exp(ll) # 1xN
    ll, norm = sess.run((ll, norm))
    return mus[:, O:D], sigmas[:, O:D], np.transpose(ll / norm)

# A "sparser" estimate which just uses the most likely cluster's mean as the estimate.
def argmax_exp(mus, sigmas, alphas):
    i = np.argmax(alphas)
    return mus[i]


# In[8]:

def MVN(shear, shift):
    rs = np.random.randn(n_samples, 2)
    return np.dot(rs, shear.T) + shift

# In[10]:

tissues_to_cluster = train_tissues
X_np = train_df_imp[train_tissues].values.transpose()

np.random.shuffle(X_np)

# In[11]:
#print(X_np.shape)

# In[17]:

# leave one out CV. Try for a given K (provided in arguments). For each fold (for a given K), perform 10 random restarts
# and pick EM fit with highest mean likelihood
# TODO: change from mean imputation

#print(X_np.shape)
X_np_new = X_np[:][:33]
X_np_new = np.delete(X_np_new, [14,25,26], 0)
#print(X_np_new.shape)

o = np.ones(len(train_df))
o[train_ix] = 0
o[test_ix] = 0
unobserved_untested_ix = np.where(o)[0]
o = np.zeros(len(train_df))
o[test_ix] = 1
unobserved_tested_ix = np.where(o)[0]

perm = np.hstack((train_ix, unobserved_tested_ix, unobserved_untested_ix))
X_np_perm = X_np_new[:, perm]
args_K = range(2, 11)

# these are global variables used in get_avg_r2_cv; don't change
#rmse_matrix = np.empty((len(args_K), X_np_perm.shape[0]))
#r2_matrix = np.empty((len(args_K), X_np_perm.shape[0]))

def get_avg_r2_cv(K):
    folds = range(0, X_np_perm.shape[0])
    r2_arr = []
    rmse_arr = []
    for fold in folds: 
        tf.reset_default_graph()
        sess = tf.InteractiveSession()
        X_test = X_np_perm[fold, :]
        X_train = np.delete(X_np_perm, fold, 0)
        #print(X_train.shape)

        lp = None
        m = None
        s = None
        a = None
        for j in range(10):
            rmu = np.random.choice(range(X_train.shape[0]), K, replace=False)
            mu_init = X_train[rmu]
            cur_lp, cur_m, cur_s, cur_a = fit_em(X_train, mu_init, 100, EPS, sess)
            if lp is None or cur_lp > lp:
                lp = cur_lp
                m = cur_m
                s = cur_s
                a = cur_a

        observed = X_test[:len(train_ix)]
        marginal_means, marginal_covs, marginal_alphas = marginal_posterior(observed.reshape(1, len(observed)), m, s, a, sess)

        pred = argmax_exp(marginal_means, marginal_covs, marginal_alphas[0])[:len(unobserved_tested_ix)]
        actual = X_test[len(train_ix):len(train_ix)+len(unobserved_tested_ix)]
        rmse = sklearn.metrics.mean_squared_error(actual, pred)
        rmse_arr.append(rmse)
        print("K: " + str(K) + ", Fold: " + str(fold))
        print('rmse:', np.sqrt(rmse)) # rmse of GMM
        r2 = sklearn.metrics.r2_score(actual, pred)
        r2_arr.append(r2)
        print('r2:', r2)
        sess.close()
    return rmse_arr, r2_arr
    #rmse_matrix[index, :] = rmse_arr
    #r2_matrix[index, :] = r2_arr

# results = Parallel(n_jobs = 2)(delayed(get_avg_r2_cv)(K) for K in args_K)
#threads = [None] * len(args_K)
#for i in range(len(args_K)):
    #threads[i] = threading.Thread(target=get_avg_r2_cv, args=(args_K[i], i))
    #threads[i].start()
#for i in range(len(args_K)):
    #threads[i].join()


print("K: " + str(K))
print("IS_MLE: " + str(IS_MLE))
print("MIN_COVAR: " + str(MIN_COVAR))
print("OUTPUT_FILE: " + OUTPUT_FILE)
rmse_arr, r2_arr = get_avg_r2_cv(K)
    
print("RMSE array")
print(rmse_arr)
print("R^2 array")
print(r2_arr)

# In[18]:

# save RMSE / R^2
np.save(OUTPUT_FILE + "_rmse", rmse_arr)
np.save(OUTPUT_FILE + "_r2", r2_arr)

