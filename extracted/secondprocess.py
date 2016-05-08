import glob, os, sys, csv
from collections import defaultdict
from itertools import groupby
import itertools
import pandas as pd
#os.chdir("/mydir")

filedf75 = pd.read_csv('GSM142175.csv', sep='\t', header = 0)
del filedf75['Unnamed: 0']
filedf76 = pd.read_csv('GSM142176.csv', sep='\t', header = 0)
del filedf76['Unnamed: 0']
filedf77 = pd.read_csv('GSM142177.csv', sep='\t', header = 0)
del filedf77['Unnamed: 0']
filedf78 = pd.read_csv('GSM142178.csv', sep='\t', header = 0)
del filedf78['Unnamed: 0']
filedf79 = pd.read_csv('GSM142179.csv', sep='\t', header = 0)
del filedf79['Unnamed: 0']
filedf80 = pd.read_csv('GSM142180.csv', sep='\t', header = 0)
del filedf80['Unnamed: 0']
filedf81 = pd.read_csv('GSM142181.csv', sep='\t', header = 0)
del filedf81['Unnamed: 0']
filedf82 = pd.read_csv('GSM142182.csv', sep='\t', header = 0)
del filedf82['Unnamed: 0']
filedf83 = pd.read_csv('GSM142183.csv', sep='\t', header = 0)
del filedf83['Unnamed: 0']
filedf84 = pd.read_csv('GSM142184.csv', sep='\t', header = 0)
del filedf84['Unnamed: 0']
filedf85 = pd.read_csv('GSM142185.csv', sep='\t', header = 0)
del filedf85['Unnamed: 0']
filedf86 = pd.read_csv('GSM142186.csv', sep='\t', header = 0)
del filedf86['Unnamed: 0']
filedf87 = pd.read_csv('GSM142187.csv', sep='\t', header = 0)
del filedf87['Unnamed: 0']
filedf88 = pd.read_csv('GSM142188.csv', sep='\t', header = 0)
del filedf88['Unnamed: 0']
filedf89 = pd.read_csv('GSM142189.csv', sep='\t', header = 0)
del filedf89['Unnamed: 0']
filedf90 = pd.read_csv('GSM142190.csv', sep='\t', header = 0)
del filedf90['Unnamed: 0']

globdf = pd.merge(filedf75, filedf76, on='start', how='inner')
print '75,76'
globdf = pd.merge(globdf, filedf77, on='start', how='inner')
print '77'
globdf.to_csv('out757677.csv', sep='\t')   
globdf = pd.merge(filedf78, filedf79, on='start', how='inner')
print '78,79'
globdf = pd.merge(globdf, filedf80, on='start', how='inner')
print '80'
globdf.to_csv('out787980.csv', sep='\t')   
globdf = pd.merge(filedf81, filedf82, on='start', how='inner')
print '81,82'
globdf = pd.merge(globdf, filedf83, on='start', how='inner')
print '83'
globdf.to_csv('out818283.csv', sep='\t')
globdf = pd.merge(filedf84, filedf85, on='start', how='inner')
print '84,85'
globdf = pd.merge(globdf, filedf86, on='start', how='inner')
print '86'
globdf.to_csv('out848586.csv', sep='\t')
globdf = pd.merge(filedf87, filedf88, on='start', how='inner')
print '87,88'
globdf = pd.merge(globdf, filedf89, on='start', how='inner')
print '89'
globdf.to_csv('out878889.csv', sep='\t')    

   


        



