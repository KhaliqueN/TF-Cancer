##############################################################################################
# Purpose: top scoring splicing events
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
prct <- 5 ## top prct percent of expression genes of perturbed splicing
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))
all_cancer <- substr(basename(all_files), 1,4)

# ## Master regulator data ---------------
# ## From this study: Predicting master transcription factors from pan-cancer expression data ----
# ## Higher CaCTS score is better ---
# tmas <- as.data.frame(readxl::read_excel('../data/sciadv.abf6123_tables_s1_to_s14.xlsx', 3))
# coln <- as.character(tmas[2,])
# wh <- which(coln %in% all_cancer)
# tmas1 <- tmas[,c(1,wh)]
# coln <- as.character(tmas1[2,])
# wh <- which(coln %in% all_cancer)
# coln <- coln[wh]
# colnames(tmas1) <- c('Symbol',coln)
# tmas2 <- tmas1[-c(1,2,3),]

# tmas <- as.data.frame(readxl::read_excel('../data/sciadv.abf6123_tables_s1_to_s14.xlsx', 6))
# tmas <- tmas[-c(1,2),]
# tmas3 <- tmas[,c(1,2,3)]
# colnames(tmas3) <- c('Cancer','Symbol','Cacts')


# tmas <- as.data.frame(readxl::read_excel('../data/sciadv.abf6123_tables_s1_to_s14.xlsx', 7))
# tmas4 <- tmas[-1,]
# colnames(tmas4) <- c('Cancer','Symbol')
# ## -----------------------------------------------------------------------

##-- Mean expression of TFs according to RNAseq data downloaded by me ---
save_dirx <- '../data/Diff_expr'
# save_dirf <- '../../../public_data/TCGA_FPKM_counts'
save_dirf <- '../../../public_data/TCGA_normalized_counts'
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
emap <- data.table::fread('../data/ensembl_name_map.txt')
emap <- unique(emap[,c(1,8)])
diffExprFiles <- list.files(save_dirx, full.names=TRUE) 
exprFiles <- list.files(save_dirf, full.names=TRUE) 
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
    # tempy_neg <- tempy[tempy$MEDIAN_DIFF < 0,]
    # tempy_pos <- tempy[tempy$MEDIAN_DIFF > 0,]
    # tempy_neg <- tempy_neg[order(tempy_neg$MEDIAN_DIFF),]
    # tempy_pos <- tempy_pos[order(-tempy_pos$MEDIAN_DIFF),]
    tempy_sorted <- tempy[order(-abs(tempy$MEDIAN_DIFF)),]
    tempy_sorted_topx <- tempy_sorted[seq(1,ceiling((prct*length(tempy[[1]]))/100)),]
    temptfs <- gtools::mixedsort(unique(tempy_sorted_topx$symbol))
   
    ##----- expression of genes ----
    # temp_filet <- data.table::fread(diffExprFiles[k])
    temp_filed <- as.data.frame(data.table::fread(exprFiles[k]))
    # temp_filed$Ensembl_gene_id <- temp_filet$Ensembl_gene_id
    temp_filex <- temp_filed[temp_filed$Ensembl_gene_id %in% tfs$Ensembl_Gene_ID, ]
    tfid <- temp_filex$Ensembl_gene_id
    last_flag <- substr(colnames(temp_filex), 14,15)
    n_col1 <- which(last_flag == '11') ## control samples
    n_col2 <- which(last_flag == '01') ## cancer samples
    temp_filex <- temp_filex[, n_col2]
    temp_filex$median <- matrixStats::rowMedians(as.matrix(temp_filex))
    temp_filex$Ensembl_gene_id <- tfid
    temp_filex_sorted <- temp_filex[order(-temp_filex$median), ]
    temp_filex_sorted <- merge(temp_filex_sorted, emap, by='Ensembl_gene_id')
    # temptfs_expr <- gtools::mixedsort(unique(temp_filex_sorted_topx$HGNC_symbol))
    temp_splice_expr <- temp_filex_sorted[temp_filex_sorted$HGNC_symbol %in% temptfs, ]

    ##------ put the median into the splice dataframe ---
    median_val <- rep(0,length(tempy_sorted_topx[[1]]))
    for(j in 1:length(temp_splice_expr[[1]])){
        wh <- which(tempy_sorted_topx$symbol == temp_splice_expr$HGNC_symbol[j])
        median_val[wh] <- temp_splice_expr$median[j]
    }
    tempy_sorted_topx$MEDIAN_EXPR <- median_val

    tcancer <- c(tcancer, rep(all_cancer[k], length(tempy_sorted_topx[[1]])))
    asid <- c(asid, tempy_sorted_topx$as_id)
    median_diff <- c(median_diff, tempy_sorted_topx$MEDIAN_DIFF)
    median_expr <- c(median_expr, tempy_sorted_topx$MEDIAN_EXPR)
    tgenes <- c(tgenes, tempy_sorted_topx$symbol)
}

##--- plot the number of splicing events affecting TFs -------------------
pdata <- data.frame(CANCER=tcancer, GENE=tgenes, ASID=asid, MEDIAN_DIFF=median_diff, MEDIAN_EXPR=median_expr)

p <- ggplot(pdata, aes(MEDIAN_DIFF, MEDIAN_EXPR, color=CANCER)) + 
geom_point()+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="Difference between median cancer\n and median normal PSI") +
scale_y_continuous(name="Median gene expression") +
# geom_text(aes(label=value), position=position_stack(vjust=0.5), size=3)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
    '#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Splice_expression.png"),width=7, height=3, dpi=400)


p <- ggplot(pdata, aes(CANCER, MEDIAN_EXPR)) + 
geom_boxplot(outlier.shape = NA)+geom_jitter(aes(color=CANCER))+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer") +
scale_y_continuous(name="Median gene expression \n in cancer samples") +
# geom_text(aes(label=value), position=position_stack(vjust=0.5), size=3)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
    '#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color='none')
ggsave(p,filename=paste0(save_dir,"/Top5_gene_expression.png"),width=7, height=3, dpi=400)

