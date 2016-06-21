import argparse
import re
import pandas as pd


def process_refseq(chrom, refseq_folder, keys):
    refs = dict(zip(keys, [[] for k in range(len(keys))]))
    with open(refseq_folder + 'human.{0}.protein.gpff'.format(chrom), 'r') as f:
        block = []
        for line in f:
            line = line.strip()
            if line != '//':
                block.append(line)
            else:
                ids = dict(zip(keys, [None for k in range(len(keys))]))
                for bl in block:
                    if bl.startswith('VERSION'):
                        ids['Protein'] = bl.split()[1].split('.')[0]
                    elif bl.startswith('DBSOURCE'):
                        ids['mRNA'] = bl.split()[-1].split('.')[0]
                    elif '/gene=' in bl:
                        ids['Gene'] = re.search('.+\"(.+)\"', bl).group(1)
                    elif '/db_xref="GeneID:' in bl:
                        ids['Entrez'] = re.search('.+\:([0-9]+)\"', bl).group(1)
                for k in keys:
                    refs[k].append(ids[k])
                block = []
    return refs


def main(args, keys):
    all_refs = dict(zip(keys, [[] for k in range(len(keys))]))
    for chrom in range(1, 27):
        refs = process_refseq(chrom, args.refseq_folder, keys)
        for k in keys:
            all_refs[k].extend(refs[k])
        print 'Done processing', args.refseq_folder + 'human.{0}.protein.gpff'.format(chrom)
    refseq = pd.DataFrame(all_refs)
    refseq.to_csv(args.output_file, sep='\t', header=True, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script creates an annotation file for RefSeq protein/transcript IDs, HUGO symbols, and Entrez gene IDs.')
    parser.add_argument('--refseq-folder', 
        metavar='<refseq-folder/>', 
        help='folder with RefSeq *.gpff files', 
        required=False)
    parser.add_argument('--output-file', 
        metavar='<output-file>', 
        help='filename for output file', 
        required=False)
    args = parser.parse_args()
    main(args, keys=['Protein', 'mRNA', 'Gene', 'Entrez'])


