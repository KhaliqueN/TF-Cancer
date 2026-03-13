##############################################################################################
# Purpose: Download human TFs and the corresponding uniprot files. 
# Download the human gene regulatory network.
# Download TCGA SpliceSeq data
##############################################################################################

rm(list=ls())

dir.create('data')

## Download transcription factor information from the paper: "The Human Transcription Factors, Cell, 2018"
xx <- "1-s2.0-S0092867418301065-mmc2.xlsx"
system(paste0("wget -O ../data/",xx," https://ars.els-cdn.com/content/image/",xx))

tf_file <- readxl::read_excel(paste0('../data/',xx),2)
tf_file <- as.data.frame(tf_file)
tf_file1 <- tf_file[,1:7]
tf_file1 <- tf_file1[-1,]
colnames(tf_file1) <- c('Ensembl_Gene_ID','Gene_Symbol','DBD','TF_indicator','TF_assessment','Binding_mode','Motif_status')

data.table::fwrite(tf_file1, 'data/filtered_TFs.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)


tf_file2 <- tf_file1[tf_file1$TF_indicator == 'Yes', ]

data.table::fwrite(tf_file2, 'data/filtered_TFs_curated.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)


# TCGA SpliceSeq data -------------------------------------

##--- Manually dowload the PSI files for these 15 cancer types from https://bioinformatics.mdanderson.org/TCGASpliceSeq/
##-- unzip the files in a specified folder --
out_folder <- 'data/PSI_data'
dir.create(out_folder)
system(paste("unzip '../data/*.zip' -d ",out_folder))
##---
