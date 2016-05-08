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
from tf_gmm_em import *

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

EPS = np.finfo(float).eps
K = int(sys.argv[1])
IS_MLE = False
if sys.argv[2] == "mle": 
  IS_MLE = True
if sys.argv[3] == "eps":
  MIN_COVAR = EPS
else:
  MIN_COVAR = float(sys.argv[3]) 
OUTPUT_FILE = sys.argv[4]


train_df_int = train_df
for i in train_tissues:
    train_df_int[i] = data_loader.local_impute(train_df[i].copy())


# In[8]:

def MVN(shear, shift):
    rs = np.random.randn(n_samples, 2)
    return np.dot(rs, shear.T) + shift

# In[10]:

tissues_to_cluster = train_tissues
X_np = train_df_int[train_tissues].values.transpose()

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
        #tf.reset_default_graph()
        #sess = tf.InteractiveSession()
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
            cur_lp, cur_m, cur_s, cur_a = fit_em(X_train, mu_init, 100, EPS, MIN_COVAR)
            if lp is None or cur_lp > lp:
                lp = cur_lp
                m = cur_m
                s = cur_s
                a = cur_a

        observed = X_test[:len(train_ix)]
        marginal_means, marginal_covs, marginal_alphas = marginal_posterior(observed.reshape(1, len(observed)), m, s, a)

        if not IS_MLE:
          pred = argmax_exp(marginal_means, marginal_covs, marginal_alphas[0])[:len(unobserved_tested_ix)]
        else:
          print('means', np.any(np.isnan(marginal_means)), np.any(np.isfinite(marginal_means)))
          print('covs', np.any(np.isnan(marginal_covs)), np.any(np.isfinite(marginal_covs)))
          print('alphas', np.any(np.isnan(np.exp(marginal_alphas[0]))), np.any(np.isfinite(np.exp(marginal_alphas[0]))))
          pred = gd_mle(marginal_means, marginal_covs, np.exp(marginal_alphas[0]), 50, EPS, warning='gd', verbose=False, minstep=1e-4)
          if pred is None:
            print('gd computation for mle screwed up, skipping')
            continue
          pred = pred.reshape(-1)[:len(unobserved_tested_ix)]
        actual = X_test[len(train_ix):len(train_ix)+len(unobserved_tested_ix)]
        rmse = sklearn.metrics.mean_squared_error(actual, pred)
        rmse_arr.append(rmse)
        print("K: " + str(K) + ", Fold: " + str(fold))
        print('rmse:', np.sqrt(rmse)) # rmse of GMM
        r2 = sklearn.metrics.r2_score(actual, pred)
        r2_arr.append(r2)
        print('r2:', r2)
        #sess.close()
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

