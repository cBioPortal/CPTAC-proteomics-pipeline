#!/bin/bash

cd $REPOS/cbioportal-test
mkdir -p download
mkdir -p data
mkdir -p refseq
mkdir -p resources

wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.*.protein.gpff.gz -P refseq/
python process_refseq.py --refseq-folder refseq/ --output-file refseq.tsv

wget -i http_breast.txt -P download/
wget -i http_ovarian.txt -P download/
wget -i http_colorectal.txt -P download/

# The following were on CDAP files used for the study
python cbio_cdap.py \
--proteome-files download/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv \
--proteome-pipeline itraq \
--ptm-files download/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv \
--ptm-pipeline itraq \
--ptm-prefix p \
--annotation refseq.tsv \
--output-file data/data_breast_protein_level.txt \
--meta-file data/meta_breast_protein_level.txt \
--cancer-id brca_tcga_pub

python cbio_cdap.py \
--proteome-files download/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv,download/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv \
--proteome-pipeline itraq \
--ptm-files download/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv \
--ptm-pipeline itraq \
--ptm-prefix p \
--annotation refseq.tsv \
--output-file data/data_ovarian_protein_level.txt \
--meta-file data/meta_ovarian_protein_level.txt \
--cancer-id ov_tcga

python cbio_cdap.py \
--proteome-files download/TCGA_Colon_VU_Proteome_CDAP.r2.precursor_area.tsv \
--proteome-pipeline precursor_area \
--annotation refseq.tsv \
--output-file data/data_colorectal_protein_level.txt \
--meta-file data/meta_colorectal_protein_level.txt \
--cancer-id coadread_tcga

# The following are on test MaxQuant files
python cbio_maxquant.py \
--proteome-files data/ibaq/proteinGroups.txt \
--ptm-files "data/ibaq/Oxidation (M)Sites.txt" \
--sample-regex std[0-9]+ \
--pipeline unlabeled \
--ptm-prefixes o \
--output-file data/maxquant_lfq_test.txt \
--meta-file data/meta_maxquant_lfq_test.txt \
--cancer-id test_id

python cbio_maxquant.py \
--proteome-files data/tmt/proteinGroups.txt \
--ptm-files "data/tmt/Oxidation (M)Sites.txt" \
--sample-regex std[0-9]+ \
--pipeline labeled \
--ptm-prefixes o \
--output-file data/maxquant_tmt_test.txt \
--meta-file data/meta_maxquant_tmt_test.txt \
--cancer-id test_id
