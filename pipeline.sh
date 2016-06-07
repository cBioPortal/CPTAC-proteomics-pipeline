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

python process_proteome.py -i download/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv -p itraq -a refseq.tsv
python process_proteome.py -i download/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv -p itraq -a refseq.tsv
python process_proteome.py -i download/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv -p itraq -a refseq.tsv
python process_proteome.py -i download/TCGA_Colon_VU_Proteome_CDAP.r2.precursor_area.tsv -p precursor_area -a refseq.tsv

python process_ptm.py -i download/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv -p itraq -a refseq.tsv
python process_ptm.py -i download/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv -p itraq -a refseq.tsv

python data_explore.py -i data/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv -t "Ovarian Proteome 1"
python data_explore.py -i data/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv -t "Ovarian Proteome 2"
python data_explore.py -i data/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv -t "Breast Proteome"
python data_explore.py -i data/TCGA_Colon_VU_Proteome_CDAP.r2.precursor_area.tsv -t "Colon Proteome"

