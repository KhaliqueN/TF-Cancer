##############################################################################################
# Purpose: Download human TFs and the corresponding uniprot files. 
# Download the human gene regulatory network.
# Download TCGA SpliceSeq data
##############################################################################################

rm(list=ls())

dir.create('../data')

## Download transcription factor information from the paper: "The Human Transcription Factors, Cell, 2018"
xx <- "1-s2.0-S0092867418301065-mmc2.xlsx"
system(paste0("wget -O ../data/",xx," https://ars.els-cdn.com/content/image/",xx))

tf_file <- readxl::read_excel(paste0('../data/',xx),2)
tf_file <- as.data.frame(tf_file)
tf_file1 <- tf_file[,1:7]
tf_file1 <- tf_file1[-1,]
colnames(tf_file1) <- c('Ensembl_Gene_ID','Gene_Symbol','DBD','TF_indicator','TF_assessment','Binding_mode','Motif_status')

data.table::fwrite(tf_file1, '../data/filtered_TFs.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)


tf_file2 <- tf_file1[tf_file1$TF_indicator == 'Yes', ]

data.table::fwrite(tf_file2, '../data/filtered_TFs_curated.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)


## Uniprot gene ID mapping file ----
uniprot_map <- data.table::fread('../data/uniprotkb_reviewed_true_AND_model_organ_2025_03_02.tsv')


## Gene regulatory network ---------------------------------
## Detected via small-scale experiments --> more confident ---
xx <- "TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv"
system(paste0("wget -O ../data/",xx," https://cdn.netbiol.org/tflink/download_files/",xx))

## Detected via large-scale (high-throughput) experiments --> relatively less confident in comparison to small-scale experiments ---
xx <- "TFLink_Homo_sapiens_interactions_LS_simpleFormat_v1.0.tsv.gz"
system(paste0("wget -O ../data/",xx," https://cdn.netbiol.org/tflink/download_files/",xx))
system(paste0("gunzip ../data/",xx))

## TCGA SpliceSeq data -------------------------------------

# ##--- Manually dowload the PSI files for these 15 cancer types from https://bioinformatics.mdanderson.org/TCGASpliceSeq/
# ##-- unzip the files in a specified folder --
# out_folder <- '../data/PSI_data'
# dir.create(out_folder)
# system(paste("unzip '../data/*.zip' -d ",out_folder))
# ##---


## Ensembl map data ------------------------------------------

##---- Download all ensembl data -------------------------
# wget -O result.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_gene_id_version" /><Attribute name = "ensembl_transcript_id" /> <Attribute name = "ensembl_transcript_id_version" /><Attribute name = "external_gene_name" /><Attribute name = "ensembl_exon_id" /><Attribute name = "external_synonym" /><Attribute name = "chromosome_name" /><Attribute name = "strand" /><Attribute name = "start_position" /><Attribute name = "end_position" /><Attribute name = "hgnc_symbol" /><Attribute name = "uniprotswissprot" /><Attribute name = "transcript_biotype" /></Dataset></Query>'
temp_file <- data.table::fread('result.txt')
colnames(temp_file) <- c('Ensembl_gene_id','Ensembl_gene_id_version','Ensembl_transcript_id','Ensembl_transcript_id_version',
    'External_gene_name','Ensembl_exon_id','External_synonym','Chromosome_name', 'Strand','CHR_start','CHR_end','HGNC_symbol',
    'Uniprotswissprot','Transcript_biotype')

data.table::fwrite(temp_file, '../data/ensembl_name_map.txt', sep='\t', row.names=FALSE, quote=FALSE)
