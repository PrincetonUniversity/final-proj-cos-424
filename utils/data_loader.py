'''Provides idempotent data loading'''

import numpy as np
import pandas as pd
import os, os.path
from itertools import *
from functools import partial
import pickle

from intervaltree import IntervalTree, Interval

chromosomes = [str(x) for x in [1, 2, 6, 7, 11]]

def chromosome_files(n):
    '''Returns list of files associated with chromosome n'''
    base = 'intersected_final_chr'
    spec = '_cutoff_20_'
    suffixes = ['train.bed', 'sample_partial.bed', 'sample_full.bed']
    return [base + n + spec + suffix for suffix in suffixes]

all_chromosomes = set(chain.from_iterable(chromosome_files(n) for n in chromosomes))
encode_file = 'wgEncodeRegTfbsClusteredV3.bed'
def have_chromosomes(): return all_chromosomes.issubset(set(os.listdir('data')))

def unzip_to_data_dir():
    '''Unzips raw data with various checks to cwd/data. Assumes bash. Idempotent.'''
    if 'methylation_imputation' not in [x for x in os.listdir('.') if os.path.isdir(x)]:
        raise Exception('Missing assignment repository in cwd')

    if not os.path.exists('data'):
        os.mkdir('data')

    if not have_chromosomes():
        os.system('cp methylation_imputation/data/*.bed.gz data/')
        os.system('gunzip data/*.gz')
        if not have_chromosomes():
            raise Exception('Error unpacking chromosomes data')

    annotations_files = {encode_file}
    def have_annotations(): return annotations_files.issubset(set(os.listdir('data')))
    if not have_annotations():
        os.system('cp methylation_imputation/annotations/*.bed.gz data/')
        os.system('gunzip data/*.gz')
        if not have_annotations():
            raise Exception('Error unpacking ENCODE data')

def read_tsv(name): return pd.read_csv(name, sep='\t', header=None)

def load_chromosome(n):
    '''Idempotent data loading. For a given chromosome n (a string).
    
    Returns (train_df, test_df, train_ix, test_ix, train_tissues, tfs)
    
    The first two are the train and test dataframes, and test_ix are the
    values in test_df['assayed'] that are missing and need to be imputed (with the
    correct answer being in test_df['filled'] in the corresponding locations.
    train_ix are the assayed (known) methylation values from limited microarray
    sampling (e.g., test_df['assayed'].iloc[train_ix] can be used for prediction of
    test_df['filled'].iloc[test_ix], and the former should be about equal to
    test_df['filled'].iloc[train_ix] (two different ways of sampling methylation).
    
    Imports genetic context and adds those columns to the parameter df, returning
    a merged one. tfs is the list of names of new transcription
    factors.
    
    train_tissues is a list of the names of columns with chromosome methylation values.
    
    Note that loading from scratch may take ~5 min (for merging genetic contexts).'''
    
    unzip_to_data_dir()
 
    train_chr = read_tsv('data/' + chromosome_files(n)[0])
    test_chr_partial = read_tsv('data/' + chromosome_files(n)[1])
    test_chr_full = read_tsv('data/' + chromosome_files(n)[2])

    unknown_chr_ix = np.where((test_chr_partial[5] == 0) & ~np.isnan(test_chr_full[4]))[0]
    known_chr_ix = np.where((test_chr_partial[5] == 1) & ~np.isnan(test_chr_partial[4]))[0]

    test_ix = unknown_chr_ix
    train_ix = known_chr_ix

    train_df = train_chr
    train_tissues = ['b' + str(x) for x in range(train_chr.shape[1] - 4)]
    train_df.columns = ['chromosome', 'start', 'end', 'strand'] + train_tissues

    test_df = test_chr_full
    test_df.columns = ['chromosome', 'start', 'end', 'strand', 'filled', 'assayed']
    test_df['missing'] = test_chr_partial[4]

    train_df['strand'] = train_df['strand'] == '+'
    test_df['strand'] = test_df['strand'] == '+'
    
    unzip_to_data_dir()
    
    # Add CRE info
    rawe = read_tsv('data/' + encode_file)
    rawe.drop([4, 5, 6, 7], axis=1, inplace=True)
    rawe.columns = ['chromosome', 'start', 'end', 'tf']
    tfs = set(rawe['tf'])

    # This may take a minute
    chr_name = 'chr' + n
    iv_chr = IntervalTree((Interval(*x[1]) for x in
                           rawe[rawe['chromosome'] == chr_name][['start', 'end', 'tf']].iterrows()))
    def inside_tf(iv, row):
        overlap = [x.data for x in iv[row['start']:row['end']]]
        return pd.Series({tf: (tf in overlap) for tf in tfs})
    # This may take 3 minutes
    train_df = train_df.merge(train_df[['start', 'end']].apply(partial(inside_tf, iv_chr), axis=1), copy=False,
                              left_index=True, right_index=True)
    
    return train_df, test_df, train_ix, test_ix, train_tissues, list(tfs)

def load_chromosome_cached(n):
    try:
        with open('data/' + n + '.p', 'rb') as r:
            t = pickle.load(r)
        return t
    except (pickle.PickleError, FileNotFoundError) as e:
        print('Encountered error while attempting to load cached result from data/' + n + '.p:')
        print('\t', str(e))
        print('Recomputing...')
        t = load_chromosome(n)
        with open('data/' + n + '.p', 'wb') as w:
            pickle.dump(t, w)
        print('Done')
        return t

def local_impute(data):
    #http://stackoverflow.com/questions/9537543/replace-nans-in-numpy-array-with-closest-non-nan-value
    mask = np.isnan(data)
    data[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])
    return data
