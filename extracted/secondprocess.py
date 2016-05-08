import glob, os, sys, csv
from collections import defaultdict
from itertools import groupby
import itertools
import pandas as pd
#os.chdir("/mydir")


global globdf

first = True

for file in glob.glob('*.csv'):
    with open(file, 'r') as csvfile:

        filedf = pd.read_csv(csvfile, sep='\t', header = 0)

        # processing one file #
        if first:
            globdf = filedf
            first = False
        else:
            globdf = pd.merge(globdf, filedf, on='start', how='inner')
        #print globdf.head()    
            

# save
globdf.to_csv('out.csv', sep='\t')        


        



