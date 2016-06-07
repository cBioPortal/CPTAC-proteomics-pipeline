import argparse
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('ggplot')

# set argparser
parser = argparse.ArgumentParser(description='Plots for data exploration of CPTAC intensity ratios.')
parser.add_argument('-i', metavar='in-file', help='/path/to/filename.ext', required=True)
parser.add_argument('-t', metavar='exp-type', help='for plot titles (i.e. "Ovarian Proteome")', required=True)
args = parser.parse_args()

# replace with args.i
data = pd.read_csv(args.i, sep='\t', header=0, index_col=0)

for col in [c for c in data.columns if 'TCGA' not in c]:
    data.drop(col, axis=1, inplace=True)

# replace with for loop
plt.figure()

for col in data.columns:
    plt.hist(data.ix[data[col] > 0, col], bins=100, histtype='stepfilled', color='#2B547E', alpha=0.2, normed=True)

plt.tick_params(axis='both', labelsize=16, labelcolor='k')
plt.ylabel('Frequency', fontsize=18, color='k')
plt.title(args.t, fontsize=18)
plt.xlim([0, 3])
plt.savefig(''.join(args.t.split()), bbox_inches='tight')



