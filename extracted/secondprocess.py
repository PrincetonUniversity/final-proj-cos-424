import glob, os, sys, csv
from collections import defaultdict
from itertools import groupby
import itertools
import pandas as pd
#os.chdir("/mydir")


global globdf
first = True
filedf75 = pd.read_csv('GSM142175.csv', sep='\t', header = 0)
filedf76 = pd.read_csv('GSM142176.csv', sep='\t', header = 0)
filedf77 = pd.read_csv('GSM142177.csv', sep='\t', header = 0)
filedf78 = pd.read_csv('GSM142178.csv', sep='\t', header = 0)
filedf79 = pd.read_csv('GSM142179.csv', sep='\t', header = 0)
filedf80 = pd.read_csv('GSM142180.csv', sep='\t', header = 0)
filedf81 = pd.read_csv('GSM142181.csv', sep='\t', header = 0)
filedf82 = pd.read_csv('GSM142182.csv', sep='\t', header = 0)
filedf83 = pd.read_csv('GSM142183.csv', sep='\t', header = 0)
filedf84 = pd.read_csv('GSM142184.csv', sep='\t', header = 0)
filedf85 = pd.read_csv('GSM142185.csv', sep='\t', header = 0)
filedf86 = pd.read_csv('GSM142186.csv', sep='\t', header = 0)
filedf87 = pd.read_csv('GSM142187.csv', sep='\t', header = 0)
filedf88 = pd.read_csv('GSM142188.csv', sep='\t', header = 0)
filedf90 = pd.read_csv('GSM142190.csv', sep='\t', header = 0)

globdf = pd.merge(filedf75, filedf76, on='start', how='inner')
print '75,76'
globdf = pd.merge(globdf, filedf77, on='start', how='inner')
print '77'
globdf = pd.merge(globdf, filedf78, on='start', how='inner')
print '79'
globdf = pd.merge(globdf, filedf79, on='start', how='inner')
print '79'
globdf = pd.merge(globdf, filedf80, on='start', how='inner')
print '80'
globdf = pd.merge(globdf, filedf81, on='start', how='inner')
print '81'
globdf = pd.merge(globdf, filedf82, on='start', how='inner')
print '82'
globdf = pd.merge(globdf, filedf83, on='start', how='inner')
print '83'
globdf = pd.merge(globdf, filedf84, on='start', how='inner')
print '84'
globdf = pd.merge(globdf, filedf85, on='start', how='inner')
print '85'
globdf = pd.merge(globdf, filedf86, on='start', how='inner')
print '86'
globdf = pd.merge(globdf, filedf87, on='start', how='inner')
print '87'
globdf = pd.merge(globdf, filedf88, on='start', how='inner')
print '88'
globdf = pd.merge(globdf, filedf89, on='start', how='inner')
print '89'
globdf = pd.merge(globdf, filedf90, on='start', how='inner')
print '90'

# save
globdf.to_csv('out.csv', sep='\t')        


        



