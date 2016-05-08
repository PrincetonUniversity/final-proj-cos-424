import glob, os, sys, csv
from collections import defaultdict
from itertools import groupby
import itertools
import pandas as pd
#os.chdir("/mydir")

global col 


def main():
    col = 0
    parse('GSM142175', col)
    col = 6
    parse('GSM142176', col)
    col = 16
    parse('GSM142177', col)
    col = 26
    parse('GSM142178', col)
    col = 36
    parse('GSM142179', col)
    col = 46
    parse('GSM142180', col)
    col = 56
    parse('GSM142181', col)
    col = 66
    parse('GSM142182', col)
    col = 76
    parse('GSM142183', col)
    col = 86
    parse('GSM142184', col)
    col = 96
    parse('GSM142185', col)
    col = 106
    parse('GSM142186', col)
    col = 116
    parse('GSM142187', col)
    col = 126
    parse('GSM142188', col)
    col = 136
    parse('GSM142189', col)
    col = 146
    parse('GSM142190', col)


def parse(prefix, col):
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

if __name__ == "__main__":
        main()
        



