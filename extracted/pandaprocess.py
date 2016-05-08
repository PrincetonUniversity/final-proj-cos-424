import glob, os, sys, csv
from collections import defaultdict
from itertools import groupby
import itertools
import pandas as pd
#os.chdir("/mydir")


global globdf

first = True
col = 0
for file in glob.glob('GSM14217*.txt'):
    print os.path.basename(file)
    with open(file, 'r') as csvfile:

        filedf = pd.read_csv(csvfile, sep='\t', header = 0)

        # processing one file #

        # take chr1
        chr1df = filedf.loc[filedf['chr'] == 'chr1']
        
        # add ratio
        ratio = ((chr1df.MethylatedReadCount)+(chr1df.ConcordantMethylatedReadCount)) / \
                ((chr1df.MethylatedReadCount)+(chr1df.UnmethylatedReadCount)+(chr1df.ConcordantMethylatedReadCount)+(chr1df.ConcordantUnmethylatedReadCount))

        chr1df['ratio' + str(col)] = ratio

        chr1df = chr1df[['start', 'ratio' + str(col)]]
        if first:
            globdf = chr1df
            first = False
        else:
            globdf = pd.merge(globdf, chr1df, on='start', how='inner')

        col = col +1
        # print col
        #print globdf.head()    
            
# save
globdf.to_csv('out7.csv', sep='\t')        


        



