##############################################################################################
# Purpose: Which TFs are expressed
##############################################################################################

rm(list=ls())
library(ggplot2)
library(data.table)

expr_cutoff <- 0.5

save_dirx <- '../data/Diff_expr'
save_dirf <- '../../../public_data/TCGA_FPKM_counts'
# save_dirf <- '../../../public_data/TCGA_normalized_counts'

tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
emap <- data.table::fread('../data/ensembl_name_map.txt')
emap <- unique(emap[,c(1,8)])


diffExprFiles <- list.files(save_dirx, full.names=TRUE) 
exprFiles <- list.files(save_dirf, full.names=TRUE) 

for(k in 1:length(diffExprFiles)){

    temp_filed <- data.table::fread(diffExprFiles[k])
    temp_filex <- temp_filed[temp_filed$Ensembl_gene_id %in% tfs$Ensembl_Gene_ID, ]
    wh1 <- which(temp_filex$log2FoldChange < -2)
    wh2 <- which(temp_filex$log2FoldChange > 2)
    temp_filey <- temp_filex[union(wh2, wh1), ]
    temp_filez <- temp_filey[temp_filey$padj < 0.05, ]

    temp_file <- as.data.frame(data.table::fread(exprFiles[k]))
    last_flag <- substr(colnames(temp_file), 14,15)
    n_col1 <- which(last_flag == '11') ## control samples
    n_col2 <- which(last_flag == '01') ## cancer samples
    temp_file <- temp_file[, n_col2]
    temp_file$mean <- rowMeans(temp_file)
    temp_file$min <- apply(temp_file, 1, FUN = min)
    temp_file$max <- apply(temp_file, 1, FUN = max)

    temp_file$Ensembl_gene_id <- temp_filed$Ensembl_gene_id
    temp_filex <- temp_file[temp_file$Ensembl_gene_id %in% temp_filez$Ensembl_gene_id, ]
    wh1 <- which(temp_filex$min > 0.1)
    temp_filey <- temp_filex[wh1, ]
    

    temp_filexx <- temp_file[temp_file$Ensembl_gene_id %in% tfs$Ensembl_Gene_ID, ]
    temp_filexx <- temp_filexx[,seq(length(temp_filexx)-3,length(temp_filexx))]
    temp_filexx <- temp_filexx[order(-temp_filexx$mean), ]

    temp_fileyy <- merge(temp_filexx, emap, by='Ensembl_gene_id')
    temp_fileyy <- temp_fileyy[order(-temp_fileyy$mean), ]

}


