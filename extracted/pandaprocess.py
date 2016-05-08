import glob, os, sys, csv
from collections import defaultdict
from itertools import groupby
import itertools
import pandas as pd
#os.chdir("/mydir")

global col 
col = 0

def main():
    parse('GSM142175')
    parse('GSM142176')
    parse('GSM142177')
    parse('GSM142178')
    parse('GSM142179')
    parse('GSM142180')
    parse('GSM142181')
    parse('GSM142182')
    parse('GSM142183')
    parse('GSM142184')
    parse('GSM142185')
    parse('GSM142186')
    parse('GSM142187')
    parse('GSM142188')
    parse('GSM142189')
    parse('GSM142190')


def parse(prefix)
    first = True
    for file in glob.glob(prefix + '*.txt'):
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
    globdf.to_csv(prefix + '.csv', sep='\t')        


        



