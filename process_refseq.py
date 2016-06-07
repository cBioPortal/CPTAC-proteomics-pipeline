import re
import pandas as pd


def process_refseq(chrom):
    refprot = []
    refmrna = []
    refgene = []
    refentz = []
    with open('refseq/human.{0}.protein.gpff'.format(chrom), 'r') as f:
        block = []
        for line in f:
            line = line.strip()
            if line != '//':
                block.append(line)
            else:
                prot_id = ''
                mrna_id = ''
                gene_nm = ''
                entz_id = ''
                for bl in block:
                    if bl.startswith('VERSION'):
                        prot_id = bl.split()[1].split('.')[0]
                    elif bl.startswith('DBSOURCE'):
                        mrna_id = bl.split()[-1].split('.')[0]
                    elif '/gene=' in bl:
                        gene_nm = re.search('.+\"(.+)\"', bl).group(1)
                    elif '/db_xref="GeneID:' in bl:
                        entz_id = re.search('.+\:([0-9]+)\"', bl).group(1)
                refprot.append(prot_id)
                refmrna.append(mrna_id)
                refgene.append(gene_nm)
                refentz.append(entz_id)
                block = []
    return refprot, refmrna, refgene, refentz


all_refprot = []
all_refmrna = []
all_refgene = []
all_refentz = []
for chrom in range(1, 27):
    refprot, refmrna, refgene, refentz = process_refseq(chrom)
    all_refprot.extend(refprot)
    all_refmrna.extend(refmrna)
    all_refgene.extend(refgene)
    all_refentz.extend(refentz)
    print 'Done', chrom

refseq = pd.DataFrame({'Protein': all_refprot, 'mRNA': all_refmrna, 'Gene': all_refgene, 'Entrez': all_refentz})
refseq.to_csv('RefSeqIDs.tsv', sep='\t', header=True, index=False)


