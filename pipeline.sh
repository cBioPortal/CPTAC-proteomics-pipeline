#!/bin/bash

cd $REPOS/cbioportal-test
mkdir -p download
mkdir -p data
mkdir -p refseq
mkdir -p resources

wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.*.protein.gpff.gz -P refseq/
python process_refseq.py -r refseq -o refseq.tsv

wget -i http_breast.txt -P download/
wget -i http_ovarian.txt -P download/
wget -i http_colorectal.txt -P download/

python cbio_proteomics.py \
--proteome-files download/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv \
--proteome-pipeline itraq \
--ptm-files download/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv \
--ptm-pipeline itraq \
--ptm-prefix p \
--annotation refseq.tsv \
--output-file data/data_protein_level.txt \
--meta-file data/meta_protein_level.txt \
--cancer-id brca_tcga_pub

