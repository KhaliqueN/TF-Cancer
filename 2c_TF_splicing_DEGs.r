##############################################################################################
# Purpose: Plots to see 
# - The dysregulation differences between genes regulated by AS affected TFs and gene regulated by non-AS affected TFs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(biomaRt)
library(GenomicDataCommons)
library(RColorBrewer)

save_dir <- '../results'

## TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
##----------------------------------------------------------------------

input_dir <- '../data/PSI_data'
fdr <- 0.05
diff <- 0.1
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI.txt', full.names=TRUE))
all_cancer <- substr(basename(all_files), 1,4)

##--- splicing events affecting TFs ----------------------------------------------------------------------
TFs_affected <- list()
TFs_naffected <- list()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')

    ## affected ---
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of TFs affected by AS
    tempy <- tempx[whx,]

    ## not affected ---
    wh1 <- which(temp$FDR > fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) == 0)
    wh <- intersect(wh1, wh2)
    tempx1 <- temp[wh, ]
    whx <- which(tempx1$symbol %in% tfs$Gene_Symbol) ## number of TFs affected by AS
    tempy1 <- tempx1[whx,]

    TFs_affected[[k]] <- unique(tempy$symbol)
    TFs_naffected[[k]] <- unique(tempy1$symbol)
}

##--- plot to compare the gene dysregulations ---------------------------------------
allfiles <- gtools::mixedsort(list.files('../data/Diff_expr', full.names=TRUE))
ensembl_gene_map <- data.table::fread('../data/ensembl_name_map.txt')

##--- grn without direction downloaded from here: https://tflink.net/download/
grn <- as.data.frame(data.table::fread('../data/TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv'))
grn_filt <- grn[grn$`Name.TF` %in% tfs$Gene_Symbol,]

grn_ls <- as.data.frame(data.table::fread('../data/TFLink_Homo_sapiens_interactions_LS_simpleFormat_v1.0.tsv'))
grn_ls_filt <- grn_ls[grn_ls$`Name.TF` %in% tfs$Gene_Symbol,]


##--- For each TF individually ----------------------------------------------------------------------
rnd_len <- 10
cut_off <- 0.05
tf_drg <- list()
prt_vals <- list()
nprt_vals <- list()

for(k in 1:length(all_cancer)){

    ##-- expression matrix --
    temp_expr <- data.table::fread(allfiles[k])
    temp_expr$Ensembl_gene_id <- unlist(lapply(strsplit(temp_expr$Ensembl_gene_id, '[.]'),'[[',1))
    atfs <- TFs_affected[[k]]

    pval_frac <- c()
    atfs_name <- c()
    prt_valsx <- list()
    nprt_valsx <- list()

    for(j in 1:length(atfs)){
        tf_prt <- grn_filt[grn_filt$`Name.TF` == atfs[j], ]
        if(nrow(tf_prt) >= 10){ ## 9 13 14 17 21.....
            tf_prt_genes <- unique(tf_prt$`Name.Target`)
            tf_prt_genes_en <- unique(ensembl_gene_map[ensembl_gene_map$HGNC_symbol %in% tf_prt_genes, ]$Ensembl_gene_id)
            prt_expr <- abs(temp_expr[temp_expr$Ensembl_gene_id %in% tf_prt_genes_en, ]$log2FoldChange)
            prt_valsx[[j]] <- prt_expr
            ##-- random expr ----
            pdiff <- c()
            rnd_exprx <- c()
            for(i in 1:rnd_len){
                rnd_expr <- abs(temp_expr[temp_expr$Ensembl_gene_id %in% sample(temp_expr$Ensembl_gene_id, length(tf_prt_genes)), ]$log2FoldChange)
                rnd_exprx <- c(rnd_exprx, rnd_expr)
                pdiff <- c(pdiff, wilcox.test(as.numeric(prt_expr), as.numeric(rnd_expr), paired=FALSE, alternative='greater')$p.value)
            }
            nprt_valsx[[j]] <- rnd_exprx
            ##-- Fraction of p-values greater than 0.05 ----------------
            pval_frac <- c(pval_frac, length(which(pdiff > cut_off))/rnd_len)
            atfs_name <- c(atfs_name, atfs[j])
        }
    }

    names(pval_frac) <- atfs_name
    tf_drg[[k]] <- pval_frac
    prt_vals[[k]] <- prt_valsx
    nprt_vals[[k]] <- nprt_valsx

    cat('Cancer',k,'of',length(all_cancer),'done\n')
}

all_sig <- lapply(tf_drg, function(x) x[x < cut_off])


###-- tile plot of sig affected TFs -------------------------------------
cancer <- c()
TF <- c()
pval <- c()
for(k in 1:length(all_sig)){
    temp <- all_sig[[k]]
    if(length(temp) > 0){
        cancer <- c(cancer, rep(all_cancer[k], length(temp)))
        TF <- c(TF, names(temp))
        pval <- c(pval, temp)
    }
}
pdata <- data.frame(Cancer=cancer, TF=TF, PVAL=pval*100)
pdata$pvalue <- cut(pdata$PVAL, breaks = c(0,0.1,0.5,1,5), include.lowest=TRUE)
cols <- brewer.pal(4,"Spectral")
cnx <- plyr::count(pdata$TF)
cnx <- cnx[order(cnx$freq), ]
cn1 <- rev(cnx$x)
p <- ggplot(pdata, aes(TF, Cancer)) + geom_tile(aes(fill=pvalue))+
  theme(legend.text=element_text(size=10))+scale_fill_manual(values=cols,drop=FALSE)
basesize <- 12
p <- p + theme_grey(base_size = basesize) + labs(y = "Cancer", x = "Transcription factor") +
  scale_x_discrete(limits = cn1) +
  scale_y_discrete()+
  guides(fill=guide_legend(title="% of random \nexperiments"))+
  theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 1,vjust=0.5, colour = "grey50"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(p,filename=paste0(save_dir,'/DEGs_TFs.png'),width=7, height=3.5, dpi=400)


# ##-- Distribution plot of the afffected TFs ------------------------------
# save_dirxx <- '../results/TF_expr'
# dir.create(save_dirxx)
# vals <- c()
# tags <- c()
# cans <- c()
# Genes <- c()
# for(k in 1:length(pdata[[1]])){

#     kk <- which(all_cancer == pdata[[1]][k])
#     temp1 <- prt_vals[[kk]]
#     temp2 <- nprt_vals[[kk]]
#     wh1 <- which(lengths(temp1) != 0)
#     wh2 <- which(tf_drg[[kk]] < cut_off)

#     for(j in 1:length(wh2)){
#         vals <- c(vals, temp1[wh1[wh2[j]]][[1]])
#         tags <- c(tags, rep('Actual', length(temp1[wh1[wh2[j]]][[1]])))
#         cans <- c(cans, rep(pdata[[1]][k], length(temp1[wh1[wh2[j]]][[1]])))
#         Genes <- c(Genes, rep(pdata[[2]][k], length(temp1[wh1[wh2[j]]][[1]])))
#         vals <- c(vals, temp2[wh1[wh2[j]]][[1]])
#         tags <- c(tags, rep('Random', length(temp2[wh1[wh2[j]]][[1]])))
#         cans <- c(cans, rep(pdata[[1]][k], length(temp2[wh1[wh2[j]]][[1]])))
#         Genes <- c(Genes, rep(pdata[[2]][k], length(temp2[wh1[wh2[j]]][[1]])))
#     }
# }

# pdatax <- data.frame(Cancer=cans, Tag=tags, Value=vals, Gene=Genes)

# for(k in 1:length(pdata[[1]])){

#     temp1 <- pdatax[pdatax$Cancer == pdata[[1]][k], ]
#     temp <- temp1[temp1$Gene == pdata[[2]][k], ]

#     if(nrow(temp) != 0){

#         p <-ggplot(temp, aes(x=Value, color=Tag)) +geom_density()
#         ggsave(p,filename=paste0(save_dirxx,'/',pdata[[1]][k],'_',pdata[[2]][k],'.png'),width=3.5, height=2.5, dpi=400)

#     }

# }









##---- over all TFs ------------------------------------------------------------------------------------
# ## ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# ## attrs <- c('ensembl_gene_id', 'uniprotswissprot', 'hgnc_symbol')
# ## 
# ## 
# pdiff <- c()
# for(k in 1:length(all_cancer)){

#     ##-- expression matrix --
#     temp_expr <- data.table::fread(allfiles[k])
#     temp_expr$FDR <- p.adjust(temp_expr$pvalue, 'fdr')
#     temp_expr$Ensembl_gene_id <- unlist(lapply(strsplit(temp_expr$Ensembl_gene_id, '[.]'),'[[',1))
#     # temp_expr <- temp_expr[temp_expr$FDR < 0.05, ]

#     tf_prt_ss <- grn_filt[grn_filt$`Name.TF` %in% TFs_affected[[k]], ] 
#     tf_prt_ss_count <- plyr::count(tf_prt_ss$`Name.Target`) ## in-degree of the genes

#     ## genes regulated by exactly one TF in the small-scale (more confident) grn set
#     wh <- which(tf_prt_ss_count$freq == 1) 
#     tf_prt_ss_genes <- tf_prt_ss_count$x[wh]

#     ## which of these genes are also found to be regulate by a TF in the large-scale (less confident) grn
#     tf_prt_ls <- grn_ls_filt[grn_ls_filt$`Name.Target` %in% tf_prt_ss_genes, ]
#     tf_prt_ls_genes <- unique(tf_prt_ls$`Name.Target`)

#     ## gene exclusively regulated by the affected TFs
#     tf_prt_genes <- setdiff(tf_prt_ss_genes, tf_prt_ls_genes)
#     tf_prt_genes_en <- unique(ensembl_gene_map[ensembl_gene_map$HGNC_symbol %in% tf_prt_genes, ]$Ensembl_gene_id)


#     tf_nprt_ss <- grn_filt[grn_filt$`Name.TF` %in% TFs_naffected[[k]], ]
#     tf_nprt_ss_count <- plyr::count(tf_nprt_ss$`Name.Target`) ## in-degree of the genes
#     # tf_nprt_genes <- unique(tf_nprt$`Name.Target`)

#     ## genes regulated by exactly one TF in the small-scale (more confident) grn set
#     wh <- which(tf_nprt_ss_count$freq == 1) 
#     tf_nprt_ss_genes <- tf_nprt_ss_count$x[wh]

#     ## which of these genes are also found to be regulate by a TF in the large-scale (less confident) grn
#     tf_nprt_ls <- grn_ls_filt[grn_ls_filt$`Name.Target` %in% tf_nprt_ss_genes, ]
#     tf_nprt_ls_genes <- unique(tf_nprt_ls$`Name.Target`)

#     ## gene exclusively regulated by the non-affected TFs
#     tf_nprt_genes <- setdiff(tf_nprt_ss_genes, tf_nprt_ls_genes)
#     tf_nprt_genes_en <- ensembl_gene_map[ensembl_gene_map$HGNC_symbol %in% tf_nprt_genes, ]$Ensembl_gene_id

#     # tf_prt_genes_en_no <- setdiff(tf_prt_genes_en, tf_nprt_genes_en)
#     # tf_nprt_genes_en_no <- setdiff(tf_nprt_genes_en, tf_prt_genes_en)

#     prt_expr <- abs(temp_expr[temp_expr$Ensembl_gene_id %in% tf_prt_genes_en, ]$log2FoldChange)
#     nprt_expr <- abs(temp_expr[temp_expr$Ensembl_gene_id %in% tf_nprt_genes_en, ]$log2FoldChange)

#     pdiff <- c(pdiff, wilcox.test(as.numeric(prt_expr), as.numeric(nprt_expr), paired=FALSE, alternative='greater')$p.value)
# }

# qdiff <- p.adjust(pdiff, 'fdr')


