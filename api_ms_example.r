#!/usr/bin/Rscript
library("cgdsr")
library("gplots")
library("RColorBrewer")
# open an instance and test
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)
# view and print all cancer studies
studies <- getCancerStudies(mycgds)
# select breast cancer study from query menu name
cancer_study <- "Breast Invasive Carcinoma (TCGA, Provisional)"
study_id <- studies[studies$name==cancer_study, "cancer_study_id"]
# view and print all patient cases
cases <- getCaseLists(mycgds, study_id)
# view and print all genetic profiles
profiles <- getGeneticProfiles(mycgds, study_id)
# select mass spec cases
case_id <- cases[cases$case_list_name=="Protein Quantification (Mass Spec)", "case_list_id"]
profile_id <- profiles[profiles$genetic_profile_name=="Protein levels (mass spectrometry by CPTAC)", "genetic_profile_id"]
# get genetic profile data from PAM50
genes <- c("UBE2T", "BIRC5", "NUF2", "CDC6", "CCNB1", "TYMS", "MYBL2", "CEP55", "MELK", "NDC80", "RRM2", "UBE2C", "CENPF", "PTTG1", "EXO1", "ORC6L", "ANLN", "CCNE1", "CDC20", "MKI67", "KIF2C", "ACTR3B", "MYC", "EGFR", "KRT5", "PHGDH", "CDH3", "MIA", "KRT17", "FOXC1", "SFRP1", "KRT14", "ESR1", "SLC39A6", "BAG1", "MAPT", "PGR", "CXXC5", "MLPH", "BCL2", "MDM2", "NAT1", "FOXA1", "BLVRA", "MMP11", "GPR160", "FGFR4", "GRB7", "TMEM45B", "ERBB2")
ms_data <- getProfileData(mycgds, genes, profile_id, case_id)
# filter out uninformative genes
ms_data[is.na(ms_data)] <- 0
# transpose and convert to matrix
ms_scaled <- as.matrix(t(ms_data))
# fetch TCGA subtypes
# URL: https://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt
subtypes <- read.table("BRCA.547.PAM50.SigClust.Subtypes.txt", sep="\t", header=TRUE, colClasses="character")
for (i in 1:nrow(subtypes)) {
    sample <- toString(subtypes[i, "Sample"])
    samp_vec <- strsplit(sample, "-", fixed=TRUE)[[1]]
    samp_end <- substr(samp_vec[4], 1, 2)
    sub_vec <- c(samp_vec[1:3], c(samp_end))
    rename <- paste(sub_vec, collapse=".")
    subtypes[i, "Sample"] <- rename
}
# map samples to subtypes, translate to colors
sample_colors <- rep("", length(colnames(ms_scaled)))
for (i in 1:ncol(ms_scaled)) {
    sample <- colnames(ms_scaled)[i]
    subtype <- subtypes[subtypes["Sample"]==sample, "PAM50"]
    if (subtype=="Basal") {
            sample_colors[i] <- "#ffd400" # yellow
    } else if (subtype=="Her2") {
            sample_colors[i] <- "#441500" # brown
    } else if (subtype=="LumA") {
            sample_colors[i] <- "#00a802" # green
    } else if (subtype=="LumB") {
            sample_colors[i] <- "#ad00aa" # purple
    }
}
# visualize in heatmap
colormap <- colorRampPalette(c("#007fff", "#007fff", "#000000", "#ff0000", "#ff0000"))(25)
png("heatmap.png", height=600, width=610)
heatmap.2(ms_scaled,
    dendrogram="both",
    distfun=function(x) dist(x, method="euclidean"),
    hclustfun=function(x) hclust(x, method="ward.D2"),
    col=colormap,
    trace="none",
    scale="none",
    labCol=rep("", nrow(ms_subset)),
    ColSideColors=sample_colors
)
dev.off()

