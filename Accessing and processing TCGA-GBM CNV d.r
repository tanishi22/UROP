#Accessing and processing TCGA-GBM CNV data 

#Installing packages
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(GenomicRanges)

#Query and download data

##ASCAT3 was selected as the workflow as it is the latest ASCAT sequencing platform and I encountered issues while downloading 
##data obtained using ASCAT2. Filtered using parameters from the GDC data portal. A total of 511 cases were obtained. 

query <- GDCquery(
    project = "TCGA-GBM",
    data.category = "Copy Number Variation", 
    data.type = "Gene Level Copy Number", 
    workflow.type = "ASCAT3", 
    legacy = FALSE)

GDCdownload(query)
cnv_data <- GDCprepare(query)

#Inspect downloaded cnv data 

##class is RangedSummarizedExperiment so Bioconductor packages like SummarizedExperiment, GRanges, and IRanges must be used. 

cnv_data 
dim(cnv_data) 
colData(cnv_data) #112 characteristics for each patient to be analysed
rowRanges(cnv_data) #GRanges object with 60632 ranges. Information provided on chromsome number, genomic ranges of each chromosome, strand, gene ID, and gene name)

#Subset data from chromosome 9. 

##The data is organised in GRanges lists, so we can specify the chromosome number and genomic range for 
##Data from three assays is available (copy_number, min_copy_number, max_copy_number). Select copy_number for analysis. s