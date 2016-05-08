import glob, os, sys, csv
from collections import defaultdict
from itertools import groupby
import itertools
import pandas as pd
#os.chdir("/mydir")



def main():
    global globdf
    global first 
    first = True
    globdf = pd.DataFrame(columns=['start'])
    parse('GSM1421', globdf, first)

    globdf.to_csv('last.csv', sep='\t')  


def parse(prefix, globdfval, firstval):
    global first 
    global globdf 
    first = firstval
    for file in sorted(glob.glob(prefix + '*.txt')):
        filename = os.path.basename(file)[:10]
        print filename
        with open(file, 'r') as csvfile:

            filedf = pd.read_csv(csvfile, sep='\t', header = 0)

            # processing one file #

            # take chr1
            chr1df = filedf.loc[(filedf['chr'] == 'chr1') & (filedf['strand'] == '+')]
            
            # add ratio
            ratio = ((chr1df.MethylatedReadCount)+(chr1df.ConcordantMethylatedReadCount)) / \
                    ((chr1df.MethylatedReadCount)+(chr1df.UnmethylatedReadCount)+(chr1df.ConcordantMethylatedReadCount)+(chr1df.ConcordantUnmethylatedReadCount))

            chr1df[filename] = ratio

            chr1df = chr1df[['start', filename]]
            if first:
                globdf = chr1df
                first = False
            else:
                globdf = pd.merge(globdf, chr1df, on='start', how='inner')
          

if __name__ == "__main__":
        main()
        



