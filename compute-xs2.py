
from matplotlib.colors import LogNorm
import tensorflow as tf
import numpy as np
import pandas as pd
from itertools import *
from joblib import Parallel, delayed
import sklearn
import math
import random
import sys
import multiprocessing
nproc = max(1, multiprocessing.cpu_count() - 1)

if 'utils' not in sys.path:
    sys.path.append('utils')

import data_loader

# Warnings

import warnings
warnings.filterwarnings('ignore')

# Idempotent, cached data retrieval script

print(data_loader.load_chromosome.__doc__)
train_df, test_df, train_ix, test_ix, train_tissues, tfs = \
    data_loader.load_chromosome_cached('1')

sess = tf.InteractiveSession()

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
print('nans in mean-imputed', np.isnan(train_df_imp[train_tissues]).sum().sum())
print('nans in interpolated', np.isnan(train_df_int[train_tissues]).sum().sum())

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
MIN_COVAR = EPS
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

# Make sure N, D are small so the numerically unstable verification code
# doesn't underflow.
N = 5
D = 10

X = np.random.normal(size=(N, D))
mu = X.mean(axis=0)
sigma = X.std(axis=0)
mus = np.array([mu, mu, mu * 2])
sigmas = np.array([sigma, sigma * 2, sigma])
K = len(sigmas)
alphas = np.random.dirichlet(np.ones(K), 1)[0]

mean_ll, resp = sess.run(estep(*(tf.constant(x) for x in (X, mus, sigmas, alphas))))
def normal_likelihoods(X, mu, sigma):
    exponent = -np.dot((X - mu[np.newaxis, :]) ** 2, 1 / sigma) / 2
    return (2 * np.pi) ** (-D / 2) * np.prod(sigma) ** (-1 / 2) * np.exp(exponent)
actual = np.array([normal_likelihoods(X, mu, sigma) for mu, sigma in zip(mus, sigmas)])
actual = sklearn.preprocessing.normalize(actual * alphas[:, np.newaxis], norm='l1', axis=0)
resp = sklearn.preprocessing.normalize(resp, norm='l1', axis=0)
print('actual likelihoods', actual)
print('log likelihoods   ', resp)
rmses = np.sqrt(sklearn.metrics.mean_squared_error(actual.T, resp.T, multioutput='raw_values'))
for i, rmse in enumerate(rmses):
    print('K={} rmse={}'.format(i, rmse))

# Similar pattern to
# https://gist.github.com/narphorium/d06b7ed234287e319f18

#todo try initializing covar to emprical cv computed from kmeans labels

# Runs up to max_steps EM iterations, stopping earlier if log likelihood improves
# less than tol.
# X should be an NxD data matrix, initial_mus should be KxD
# max_steps should be an int, tol should be a float.
def fit_em(X, initial_mus, max_steps, tol):
    N, D = X.shape
    K, Dmu = initial_mus.shape
    assert D == Dmu
        
    mus0 = initial_mus
    sigmas0 = np.tile(np.var(X, axis=0), (K, 1))
    alphas0 = np.ones(K)
    X = tf.constant(X)
    
    mus, sigmas, alphas = (tf.Variable(x, dtype='float64') for x in [mus0, sigmas0, alphas0])
        
    all_ll, resp = estep(X, mus, sigmas, alphas)
    cmus, csigmas, calphas = mstep(X, resp)
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
            print('EM iteration', i, 'log likelihood', ll)
            if abs(ll - prev_ll) < tol:
                print('Converged in', i, 'iterations')
                break
            prev_ll = ll
        m, s, a = sess.run((mus, sigmas, alphas))
    
    return ll, m, s, a

n_samples = 25
np.random.seed(0)

def MVN(shear, shift):
    rs = np.random.randn(n_samples, 2)
    return np.dot(rs, shear.T) + shift
G1 = MVN(np.identity(2), np.array([20, 20]))
G2 = MVN(np.array([[0., 3.5], [-0.7, .7]]), np.zeros(2))
G3 = MVN(np.array([[-1, 0.8], [0, 4]]), np.array([10, 8]))

# concatenate the two datasets into the final training set
X = np.vstack([G1, G2, G3])
rx = np.random.choice(range(len(X)), 3, replace=False)

# http://scikit-learn.org/stable/auto_examples/mixture/plot_gmm_pdf.html
_, m, s, a = fit_em(X, X[rx], 100, EPS)
x = np.linspace(-20.0, 30.0)
y = np.linspace(-20.0, 40.0)
xx, yy = np.meshgrid(x, y)
pts = np.array([xx.ravel(), yy.ravel()]).T
ll = -estep(pts, m, s, a)[0].eval()

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
    D = len(mu)
    observed_mus, observed_sigmas = (tf.constant(a, dtype='float64')
                                     for a in (mus[:,0:O], sigmas[:, 0:O]))
    observed_likelihood = tf.exp(tf_log_likelihood(xs, observed_mus, observed_sigmas, alphas)).eval() # KxN
    sklearn.preprocessing.normalize(observed_likelihood, norm='l1', axis=0, copy=False)
    return mus[:, O:D], sigmas[:, O:D], observed_likelihood.T

# http://stackoverflow.com/questions/33702251/tensorflow-loss-minimization-type-error
class DoubleGDOptimizer(tf.train.GradientDescentOptimizer):
    def _valid_dtypes(self):
        return set([tf.float32, tf.float64])

# Given a GMM defined by KxD mus, KxD sigmas, and K alphas, returns the MLE estimate
# by using gradient descent starting from each mu in mus
def gd_mle(mus, sigmas, alphas, nsteps, tol, warning=None, verbose=False, sess=None, minstep=1e-1):
    # Use GD from each of the means with a step size less than the min distance between those means
    # step size is min(mean diff) * min(mean likelihoods) / sum(adj likelihoods)
    mudist = sklearn.metrics.pairwise.pairwise_distances(mus, metric='l2')
    alphasM = alphas.reshape(-1, 1)
    alphadist = sklearn.metrics.pairwise.pairwise_distances(alphasM, metric='l1')
    min_alpha = np.minimum(np.tile(alphasM, len(alphas)), alphas)
    steps_per_decay = max(nsteps // 100, 1)
    step_size = mudist * min_alpha / alphadist / steps_per_decay
    np.fill_diagonal(step_size, np.inf)
    step_size = step_size.min(axis=1)
    best_ll = np.inf, None
    
    at_least_one_conv_early = False
    
    for i, (mu, step) in enumerate(zip(mus, step_size)):
        if verbose: print('Running mean', i)
        decay_fraction, decay_period = 0.95, steps_per_decay
        global_step = tf.Variable(0, trainable=False)
        step = max(step, minstep)
        learning_rate = tf.train.exponential_decay(
            step, global_step, decay_period, decay_fraction, staircase=True)
        learning_rate = tf.cast(learning_rate, 'float64')
        x = tf.Variable(mu.reshape(1, -1), dtype='float64')
        nlog_likelihood = -estep(x, mus, sigmas, alphas)[0]
        nlog_likelihood = tf.squeeze(nlog_likelihood)
        train_step = DoubleGDOptimizer(learning_rate).minimize(
            nlog_likelihood, global_step=global_step, var_list=[x])
        
        def get_best(sess):
            prev_ll = -np.inf
            ll = prev_ll
            sess.run(tf.initialize_all_variables())
            for i in range(1, 1 + nsteps):
                sess.run(train_step)
                ll = sess.run(nlog_likelihood)
                if verbose and i % max(nsteps // 10, 1) == 0:
                    print('  ', i, 'of', nsteps, 'done, NLL =', ll)
                if abs(ll - prev_ll) < tol:
                    if verbose: print('  Converged early after {} steps.'.format(i))
                    break
                prev_ll = ll
            if verbose: print('  Best NLL =', ll)
            return (ll, sess.run(x)), i < 1 + nsteps
        
        if not sess:
            with tf.Session() as sess:
                other, early = get_best(sess)
                at_least_one_conv_early |= early
                best_ll = min(best_ll, other)
            sess = None
        else: best_ll = min(best_ll, get_best(sess))
            
    if not at_least_one_conv_early and warning:
        print('Warning, none converged early:', warning)
        
    return best_ll[1]

# A "sparser" estimate which just uses the most likely cluster's mean as the estimate.
def argmax_exp(mus, sigmas, alphas):
    i = np.argmax(alphas)
    return mus[i]

# Make this into a nice image.
ys = [-10, 0, 2.5, 10, 20]
pts = np.arange(-20, 30, 0.1, dtype='float64').reshape(-1, 1)
ys = np.array(ys, dtype='float64').reshape(-1, 1)
# Reverse m, s, a, since we know 'y'
mr, sr = (x[:, ::-1] for x in (m, s))
mm, sm, ams = marginal_posterior(ys, mr, sr, a)
for i, y in enumerate(ys.reshape(-1)):
    print('y =', y)
    ll = estep(pts, mm, sm, ams[i])[0].eval().reshape(-1)
    #plt.show()

ys = np.arange(-20, 30, 0.05, dtype='float64').reshape(-1, 1)
mm, sm, ams = marginal_posterior(ys, mr, sr, a)
xs = [argmax_exp(mm, sm, am) for am in ams]


def f(i, am):
    ret = gd_mle(mm, sm, am, 500, EPS, warning=str(i), verbose=True)
    print('done with {}\n'.format(i+1), end='')
    sys.stdout.flush()
    return ret[0][0]


lo = int(sys.argv[1])
hi = int(sys.argv[2])

xs2 = [f(i, am) for i, am in enumerate(ams[lo:hi])]

print(xs2)

import pickle

with open('xs2-{}-{}.p'.format(lo, hi), 'wb') as w:
    pickle.dump(xs2, w)
