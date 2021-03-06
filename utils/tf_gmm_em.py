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

import warnings
warnings.filterwarnings('ignore')

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
    logc =  \
        - tf.reduce_sum(invsig * tf.square(mus), 1, keep_dims=True) \
        + 2 * tf.matmul(invsig * mus, XT) \
        - tf.matmul(invsig, tf.square(XT)) # KxN
    
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
MIN_COVAR_DEFAULT = EPS

# X is NxD, resp is KxN (and normalized along axis 0)
# Returns maximize step means, covariance, and cluster priors,
# which have dimension KxD, KxD, and K, respectively
def mstep(X, resp, min_covar=MIN_COVAR_DEFAULT):
    weights = tf.reduce_sum(resp, 1) # K
    invweights = tf.expand_dims(tf.inv(weights + 10 * EPS), 1) # Kx1
    alphas = EPS + weights / (tf.reduce_sum(weights) + 10 * EPS) # K
    weighted_cluster_sum = tf.matmul(resp, X) # KxD 
    mus = weighted_cluster_sum * invweights
    avg_X2 = tf.matmul(resp, tf.square(X)) * invweights
    avg_mu2 = tf.square(mus)
    avg_X_mu = mus * weighted_cluster_sum * invweights
    sigmas = avg_X2 - 2 * avg_X_mu + avg_mu2 + min_covar
    # (x - mu) (x-mu)^T for banded. 
    return mus, sigmas, alphas

# Similar pattern to
# https://gist.github.com/narphorium/d06b7ed234287e319f18

# Runs up to max_steps EM iterations, stopping earlier if log likelihood improves
# less than tol.
# X should be an NxD data matrix, initial_mus should be KxD
# max_steps should be an int, tol should be a float.
def fit_em(X, initial_mus, max_steps, tol, min_covar=MIN_COVAR_DEFAULT):
    tf.reset_default_graph()
    
    N, D = X.shape
    K, Dmu = initial_mus.shape
    assert D == Dmu
        
    mus0 = initial_mus
    sigmas0 = np.tile(np.var(X, axis=0), (K, 1))
    alphas0 = np.ones(K) / K
    X = tf.constant(X)
    
    mus, sigmas, alphas = (tf.Variable(x, dtype='float64') for x in [mus0, sigmas0, alphas0])
    
    all_ll, resp = estep(X, mus, sigmas, alphas)
    cmus, csigmas, calphas = mstep(X, resp, min_covar=min_covar)
    update_mus_step = tf.assign(mus, cmus)
    update_sigmas_step = tf.assign(sigmas, csigmas)
    update_alphas_step = tf.assign(alphas, calphas)     
    
    init_op = tf.initialize_all_variables()
    ll = prev_ll = -np.inf

    with tf.Session() as sess:
        sess.run(init_op)
        for i in range(max_steps):
            ll = sess.run(tf.reduce_mean(all_ll))
            sess.run((update_mus_step, update_sigmas_step, update_alphas_step))
            #print('EM iteration', i, 'log likelihood', ll)
            if abs(ll - prev_ll) < tol:
                break
            prev_ll = ll
        m, s, a = sess.run((mus, sigmas, alphas))
    
    return ll, m, s, a

# Given a set of partial observations xs each of dimension O < D for a fitted GMM model with 
# K cluster priors alpha, KxD means mus, and KxD diagonal covariances sigmas,
# returns the weighted sum of normals for the remaining D - O dimensions.
#
# Returns posterior_mus, posterior_sigmas, posterior_prior,
# of dimensions:
# Kx(D-O), Kx(D-O), NxK, respectively (each mu, sigma is the same for all posteriors).
# NxK, NxKxD, NxKxD, respectively, for each x in xs, total of N.
def marginal_posterior(xs, mus, sigmas, alphas):
    # https://gbhqed.wordpress.com/2010/02/21/conditional-and-marginal-distributions-of-a-multivariate-gaussian/
    # diagonal case is easy:
    # https://en.wikipedia.org/wiki/Schur_complement#Applications_to_probability_theory_and_statistics
    O = xs.shape[1]
    D = mus.shape[1]
    observed_mus, observed_sigmas = (tf.constant(a, dtype='float64')
                                     for a in (mus[:,0:O], sigmas[:, 0:O]))
    ll = tf_log_likelihood(xs, observed_mus, observed_sigmas, alphas) # KxN
    norm = tf_log_sum_exp(ll) # 1xN
    with tf.Session() as sess:
        ll, norm = sess.run((ll, norm))
    return mus[:, O:D], sigmas[:, O:D], np.transpose(ll - norm)

# A "sparser" estimate which just uses the most likely cluster's mean as the estimate.
def argmax_exp(mus, sigmas, alphas):
    i = np.argmax(alphas)
    return mus[i]


# http://stackoverflow.com/questions/33702251/tensorflow-loss-minimization-type-error
class DoubleGDOptimizer(tf.train.GradientDescentOptimizer):
    def _valid_dtypes(self):
        return set([tf.float32, tf.float64])

# Given a GMM defined by KxD mus, KxD sigmas, and K alphas, returns the MLE estimate
# by using gradient descent starting from each mu in mus
def gd_mle(mus, sigmas, alphas, nsteps, tol, warning=None, verbose=False, minstep=1e-3):
    # Use GD from each of the means with a step size less than the min distance between those means
    # step size is min(mean diff) * min(mean likelihoods) / sum(adj likelihoods)
    mudist = sklearn.metrics.pairwise.pairwise_distances(mus, metric='l2')
    alphasM = np.tile(alphas, (len(alphas), 1)).transpose()
    alpha_sum = np.array([[i+j for j in alphas] for i in alphas])
    steps_per_decay = max(nsteps // 100, 1)
    step_size = mudist * alphasM / alpha_sum / max(steps_per_decay, 2)
    np.fill_diagonal(step_size, np.inf)
    step_size = step_size.min(axis=1)
    best_ll = np.inf, None
    
    at_least_one_conv_early = False
    
    for i, (mu, step) in enumerate(zip(mus, step_size)):
        if verbose: print('Running mean', i)

        tf.reset_default_graph()
        decay_fraction, decay_period = 0.95, steps_per_decay
        global_step = tf.Variable(0, trainable=False, name='global_step')
        if np.isnan(step) or minstep > step: step = minstep
                    
        learning_rate = tf.train.exponential_decay(
            step, global_step, decay_period, decay_fraction, staircase=True)
        learning_rate = tf.cast(learning_rate, 'float64')
        x = tf.Variable(mu.reshape(1, -1), dtype='float64', name='x')
        nlog_likelihood = -estep(x, mus, sigmas, alphas)[0]
        nlog_likelihood = tf.squeeze(nlog_likelihood)
        train_step = DoubleGDOptimizer(learning_rate).minimize(
            nlog_likelihood, global_step=global_step, var_list=[x])
        
        with tf.Session() as session:
            prev_ll = np.inf
            ll = prev_ll
            session.run(tf.initialize_all_variables())
            for i in range(1, 1 + nsteps):
                session.run(train_step)
                ll = session.run(nlog_likelihood)
                if verbose and i % max(nsteps // 10, 1) == 0:
                    print('  ', i, 'of', nsteps, 'done, NLL =', ll)
                if abs(ll - prev_ll) < tol:
                    if verbose: print('  Converged early after {} steps.'.format(i))
                    break
                prev_ll = ll
            if verbose: print('  Best NLL =', ll)
            at_least_one_conv_early |= i < 1 + nsteps
            best_ll = min(best_ll, (ll, session.run(x)))
            
    if not at_least_one_conv_early and warning:
        print('Warning, none converged early:', warning)
        
    return best_ll[1].reshape(-1) if best_ll[1] is not None else None
