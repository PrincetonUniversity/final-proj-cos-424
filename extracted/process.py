import glob, os, sys, csv
from collections import defaultdict
from itertools import groupby
import itertools
import pandas as pd
#os.chdir("/mydir")


global globdf
globdf = pd.DataFrame(columns=['start','strand','ratio'])

for file in glob.glob('*.txt'):
    with open(file, 'r') as csvfile:
        columns =  defaultdict(list) #lists all columns and saves them in dict
        all_file = pd.read_csv(csvfile, sep='\t', header = 0)
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            for (k,v) in row.items(): 
                columns[k].append(v)

    # processing one file #

    # take chr1    
    chrlist = [len(list(group)) for key, group in groupby(columns['chr'])]
    idx = chrlist[0]
    for key, value in columns.items():
        columns[key] = value[:idx] 

    # calculate methylation values
    m_ratio = []

    for mc, uc, cmc, cuc in itertools.izip(columns['MethylatedReadCount'], columns['UnmethylatedReadCount'], columns['ConcordantMethylatedReadCount'], columns['ConcordantUnmethylatedReadCount']):
        m_ratio.append((float(mc)+float(cmc))/(float(mc)+float(uc)+float(cmc)+float(cuc)))

    subcol = defaultdict(list)

    subcol['start'] = columns['start']
    subcol['strand'] = columns['strand']
    subcol['ratio'] = m_ratio

    df = pd.DataFrame(subcol)

    pd.merge(globdf, df, on='start', how='inner')

print globdf.head()    


#output
'''
with open('out.csv', 'wb') as csvfile:
    outwriter = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
    for row in range(len(idxdict['start'])):
        outwriter.writerow([idxdict['start'][row], idxdict['strand'][row], [sublist[row] for sublist in globratio]])
'''

