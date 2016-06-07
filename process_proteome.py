import argparse
import re
import numpy as np
import pandas as pd

# set argparser
parser = argparse.ArgumentParser(description='This script converts CPTAC intensities data to importable text files for cBioPortal.')
parser.add_argument('-i', metavar='in-file', help='/path/to/filename.ext', required=True)
parser.add_argument('-p', metavar='pipeline', choices=('itraq', 'precursor_area'), help='processing pipeline (i.e. "itraq" for Breast and Ovarian, "precursor_area for Colon")', required=True)
args = parser.parse_args()

exp_fname = args.i
sample_regex = '([A-Z0-9]{2}\-[A-Z0-9]{4}\-[A-Z0-9]{2})[A-Z0-9\-]+'

# load data
exp = pd.read_csv(exp_fname, sep='\t', header=0, index_col=0)

for col in ['Mean', 'Median', 'StdDev']:
    if col in exp.index:
        exp.drop([col], inplace=True)

unique_samples = list(set([re.match(sample_regex, col).group(1) for col in exp.columns if re.match(sample_regex, col)]))

# handle different pipelines
if args.p == 'itraq':
    new_exp = np.exp2(exp[[col for col in exp.columns if re.match(sample_regex + ' Log Ratio', col)]])
elif args.p == 'precursor_area':
    new_cols = [col for col in exp.columns if re.match(sample_regex + ' Area', col)]
    new_exp = exp[new_cols].div(exp[new_cols].mean(axis=1), axis=0)

# annotate phosphosites with hugo
prot_anno = pd.read_csv('RefSeqIDs.tsv', sep='\t', header=0, dtype='str')
prot_anno = prot_anno.drop_duplicates('Gene')
gene2entz = dict(zip(prot_anno['Gene'], prot_anno['Entrez']))

entz = []
for item in new_exp.index:
    try:
        entz.append(gene2entz[item])
    except KeyError:
        entz.append('NA')

# format and save
new_exp.columns = ['TCGA-'+re.match(sample_regex, col).group(1) for col in new_exp.columns]
new_exp.index.names = ['Hugo_Symbol']
new_exp.insert(loc=0, column='Entrez_Gene_Id', value=np.array(entz))
new_exp.fillna(0, inplace=True)
new_exp.to_csv(exp_fname.replace('download/', 'data/'), sep='\t', header=True, index=True)

