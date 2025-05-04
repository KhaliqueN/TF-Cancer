##############################################################################################
# Purpose: perturbed events whose genes are NOT DEGs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)

save_dir <- '../results'

##-- TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
##----------------------------------------------------------------------

input_dir <- '../data/PSI_data'
fdr <- 0.05
diff <- 0.1
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))
all_cancer <- substr(basename(all_files), 1,4)

##-- Mean expression of TFs according to RNAseq data downloaded by me ---
save_dirx <- '../data/Diff_expr'
# save_dirf <- '../../../public_data/TCGA_FPKM_counts'
# save_dirf <- '../../../public_data/TCGA_normalized_counts'
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
emap <- data.table::fread('../data/ensembl_name_map.txt')
emap <- unique(emap[,c(1,8)])
diffExprFiles <- list.files(save_dirx, full.names=TRUE) 
# exprFiles <- list.files(save_dirf, full.names=TRUE) 
##------------------------------------------------------------------------


##--- overlap with TFs showing significant splicing events ---------------asid <- c()
median_diff <- c()
median_expr <- c()
tcancer <- c()
tgenes <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    
    ##----- diff expression of genes ----
    temp_filed <- as.data.frame(data.table::fread(diffExprFiles[k]))
    temp_filed <- temp_filed[temp_filed$pvalue > 0.05, ]
    temp_filex <- temp_filed[temp_filed$Ensembl_gene_id %in% tfs$Ensembl_Gene_ID, ]
    temp_filey <- merge(temp_filex, emap, by='Ensembl_gene_id')
    temp_filey <- temp_filey[temp_filey$HGNC_symbol != '', ]
    temp_splice_expr <- temp_filey[temp_filey$HGNC_symbol %in% tempy$symbol, ]

    ##------ put the diff expression pvalue into the splice dataframe ---
    log2f <- rep(NA,length(tempy[[1]]))
    expr_pval <- rep(NA,length(tempy[[1]]))

    for(j in 1:length(temp_splice_expr[[1]])){
        wh <- which(tempy$symbol == temp_splice_expr$HGNC_symbol[j])
        log2f[wh] <- temp_splice_expr$log2FoldChange[j]
        expr_pval[wh] <- temp_splice_expr$pvalue[j]
    }

    tempy$DEG_LOG2 <- log2f
    tempy$DEG_PVALUE <- expr_pval
    wh <- which(!is.na(tempy$DEG_LOG2))
    tempz <- tempy[wh, ]

    tempz <- tempz[order(-abs(tempz$MEDIAN_DIFF)), ]


    # tcancer <- c(tcancer, rep(all_cancer[k], length(tempy_sorted_topx[[1]])))
    # asid <- c(asid, tempy_sorted_topx$as_id)
    # median_diff <- c(median_diff, tempy_sorted_topx$MEDIAN_DIFF)
    # median_expr <- c(median_expr, tempy_sorted_topx$MEDIAN_EXPR)
    # tgenes <- c(tgenes, tempy_sorted_topx$symbol)
}

# ##--- plot the number of splicing events affecting TFs -------------------
# pdata <- data.frame(CANCER=tcancer, GENE=tgenes, ASID=asid, MEDIAN_DIFF=median_diff, MEDIAN_EXPR=median_expr)

# p <- ggplot(pdata, aes(MEDIAN_DIFF, MEDIAN_EXPR, color=CANCER)) + 
# geom_point()+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_continuous(name="Difference between median cancer\n and median normal PSI") +
# scale_y_continuous(name="Median gene expression") +
# # geom_text(aes(label=value), position=position_stack(vjust=0.5), size=3)+
# scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
#     '#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(color=guide_legend(title="Cancer type",ncol=2))
# ggsave(p,filename=paste0(save_dir,"/Splice_expression.png"),width=7, height=3, dpi=400)

