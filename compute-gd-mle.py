'''
Usage: python compute-gd-mle.py full-models.p outfile lo hi

Assuming that full-models.p is a pickle file with a tuple
(mus, sigmas, alphass) for KxD means mus, KxD diagonal covariances sigmas,
and NxK cluster posteriors alphass (NOT log posteriors).

Given an interval [lo, hi), a subset of [1, N), computes pickles to outfile
the (hi-lo)xD MLE estimates for each alphass in the interval [lo, hi)
'''
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

import warnings
warnings.filterwarnings('ignore')

if len(sys.argv) != 5:
    print(__doc__)
    sys.exit(1)

import pickle

with open(sys.argv[1], 'rb') as r:
    mus, sigmas, alphass = pickle.load(r)

lo = int(sys.argv[3])
hi = int(sys.argv[4])

mles = [gd_mle(mus, sigmas, alphass[i], 100, 1e-3, warning=str(i), verbose=True,
               minstep=1e-3) for i in range(lo, hi)]

with open(sys.argv[2], 'wb') as w:
    pickle.dump(mles, w)
