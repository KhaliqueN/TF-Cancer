##############################################################################################
# Purpose: Correlation between patient-specific delta PSI values of TFs and patient-specific 
# gene expression values of individual genes
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)
library(dorothea)
library(pheatmap)
library(igraph)
library(Hmisc)

result_dir <- '../results'
psi_input <- '../data/PSI_data'
store_dir <- '../data/correlations'
# if(dir.exists(store_dir)){
#     unlink(store_dir, recursive=TRUE)
# }
# dir.create(store_dir)

fdr <- 0.05

rnd_expr <- 1000
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
tf_ensemb_map <- as.data.frame(data.table::fread('../data/TF_ensembl_uniprot.txt', sep='\t'))
uniprot_ensembl_map <- data.table::fread('../data/ensembl_name_map.txt')
symbol_ensembl_map <- unique(uniprot_ensembl_map[,c(1,8)])
symbol_ensembl_map <- symbol_ensembl_map[which(symbol_ensembl_map$HGNC_symbol != ''),]

all_filesxx <- gtools::mixedsort(list.files(psi_input, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_psi <- list.files('../data/PSI_data', pattern='^PSI_download', full.names=TRUE)
all_cancer <- substr(basename(all_filesxx), 1,4)


all_files_genes <- gtools::mixedsort(list.files('../../../public_data/TCGA_normalized_counts', full.names=TRUE))
## Some of the event coordinates from TCGA splice seq does not overlap to any of the exons mapped from the canonical protein
## sequence from Uniprot. This is because the canonical seqeunce is not always the longest.

##--- grn without direction downloaded from here: https://tflink.net/download/
grn <- as.data.frame(data.table::fread('../data/TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv'))
grn_filt <- grn[grn$`Name.TF` %in% tfs$Gene_Symbol,]
##--- 14,439 edges ----
genes2consider <- union(grn_filt$`Name.TF`, grn_filt$`Name.Target`)

for(k in 4:length(all_cancer)){

    temp <- data.table::fread(all_filesxx[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    if(length(wh) == 0){next}
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    txg <- tempx[whx,]$as_id

    tempxx <- as.data.frame(data.table::fread(all_psi[k], sep='\t', fill=TRUE))
    tempy <- tempxx[tempxx$as_id %in% txg,]

    wh <- which(colnames(tempy) %like% '_Norm')
    tempy1 <- tempy[, -wh]
    wh <- which(colnames(tempy1) %like% 'TCGA')

    tempy2 <- tempy1[, wh]

    ##--- gene expression ----
    temp <- as.data.frame(data.table::fread(all_files_genes[k], sep='\t'))
    labs <- temp[length(temp)]
    temp1 <- cbind(labs,temp[,-length(temp)])
    tgenes <- unlist(lapply(strsplit(temp1[[1]],'[.]'),'[[',1))
    ## keep only the Ensembl ids present in the Ensembl gene symbol mapping downloaded from Ensembl
    whg <- which(tgenes %in% symbol_ensembl_map$Ensembl_gene_id)
    temp <- temp1[,-1]
    temp <- temp[whg,]

    ttgenes <- symbol_ensembl_map[which(symbol_ensembl_map$Ensembl_gene_id %in% tgenes[whg]),]
    ttgenes <- ttgenes[match(tgenes[whg], ttgenes$Ensembl_gene_id),]
    ## filter to only keep gene in the GRN --
    whgg <- which(ttgenes$HGNC_symbol %in% genes2consider)
    temp <- temp[whgg, ]
    ttgenesx <- ttgenes[whgg,]

    flag <- substr(colnames(temp),14,15)
    why <- which(flag == '11') 
    normal_temp <- temp[,why]
    normal_temp$MEDIAN <- matrixStats::rowMedians(as.matrix(normal_temp))
    cancer_temp <- temp[,setdiff(seq(1,length(temp)), why)]
    colnames(cancer_temp) <- substr(colnames(cancer_temp), 1,12)
    colnames(cancer_temp) <- gsub('-','_',colnames(cancer_temp))

    ## compute gene counts values with respect to the median normal gene counts
    for(j in 1:length(cancer_temp)){
        txc <- cancer_temp[[j]]
        if(length(which(!is.na(as.numeric(txc)) == FALSE) > 0)){
            break
        }
        txc <- as.numeric(txc) - normal_temp$MEDIAN
        cancer_temp[j] <- txc
    }

    ## only keep the same samples in TCGA gene expression and TCGA spliceseq data
    keep <- intersect(colnames(tempy2), colnames(cancer_temp))
    tempy3 <- tempy2[, which(colnames(tempy2) %in% keep)]
    cancer_temp_filt <- cancer_temp[, which(colnames(cancer_temp) %in% keep)]

    tempy3 <- tempy3[, order(names(tempy3))]
    cancer_temp_filt <- cancer_temp_filt[, order(names(cancer_temp_filt))]
    

    ##-- check if the columns in both PSI and expression are same and in same order
    xx <- all(colnames(tempy3) %in% colnames(cancer_temp_filt)) ## all row and colnames present
    yy <- all(colnames(tempy3) == colnames(cancer_temp_filt)) ## all row and colnames present in the same order

    if(!xx | !yy){
        print('Column names of the two dataframes either not same or not in the same order or both')
        break
    }

    ##--- correlation between the psi values and the gene expressions -----
    tempy3x <- as.data.frame(sapply(tempy3,as.numeric))
    tcount <- apply(tempy3, 1, function(x) length(which(!is.na(as.numeric(x))))) 

    temp_cor <- c()
    temp_pval <- c()
    temp_rnd <- c()
    temp_as <- c()
    temp_gene <- c()
    for(j in 1:length(tempy1[[1]])){
        temp_grn <- grn_filt[grn_filt$`Name.TF` == tempy1$symbol[j], ]
        if(nrow(temp_grn) != 0){
            txc <- tempy3x[j,]
            wh <- which(!is.na(txc))
            txc <- as.numeric(txc[wh])
            for(i in 1:nrow(temp_grn)){
                wht <- which(ttgenesx$HGNC_symbol == temp_grn$`Name.Target`[i])
                if(length(wht) != 0){
                    gxc <- cancer_temp_filt[wht,]
                    gxc <- as.numeric(gxc[wh])
                    tcor <- cor.test(x=gxc, y=txc, method = 'spearman')    
                    temp_cor <- c(temp_cor, tcor$estimate)
                    temp_as <- c(temp_as, tempy1$as_id[j])
                    temp_gene <- c(temp_gene, temp_grn$`Name.Target`[i])
                    temp_pval <- c(temp_pval, tcor$p.value)
                    t_rnd_cors <- c()
                    for(i in 1:rnd_expr){ ##--- randomized experiment
                        temp_cor_rnd <- cor.test(x=gxc, y=sample(txc), method = 'spearman')
                        t_rnd_cors <- c(t_rnd_cors, abs(as.numeric(temp_cor_rnd$estimate)))
                    }
                    temp_rnd <- c(temp_rnd, length(which(t_rnd_cors > abs(as.numeric(tcor$estimate)))))
                }
            }
        }else{
            temp_cor <- c(temp_cor, NA)
            temp_as <- c(temp_as, tempy1$as_id[j])
            temp_gene <- c(temp_gene, NA)
            temp_pval <- c(temp_pval, NA)
            temp_rnd <- c(temp_rnd, NA)
        }
    }

    tempd <- data.frame(ASID=temp_as, GENE=temp_gene, COR=temp_cor, CPVAL=temp_pval, RPVAL=temp_rnd/rnd_expr)
    tempd <- tempd[complete.cases(tempd),]
    tempd$CFDR <- p.adjust(tempd$CPVAL, 'fdr')
    tempd$RFDR <- p.adjust(tempd$RPVAL, 'fdr')

    data.table::fwrite(tempd, paste0(store_dir,'/',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)
    cat('Cancer',k, 'of',length(all_cancer),'done\n')

}

# ##-- to be run ---
# # ### *****  Number of event-gene pairs **********************************************************************
# # # *****************************************************************************
# # # *****************************************************************
# # # ***************************************************

# ##--------------------
# all_sigs <- gtools::mixedsort(list.files(store_dir, pattern='*.txt', full.names=TRUE))
# num_sig_cors <- c()
# num_sig_cors_d <- list()
# background_as <- igraph::make_empty_graph(n = 0, directed = FALSE)

# for(k in 1:length(all_cancer)){
#     temp <- data.table::fread(all_sigs[k], sep='\t')
#     temp1 <- temp[temp$RFDR < 0.05, ]
#     temp2 <- temp1[temp1$CFDR < 0.05, ]
#     temp3 <- temp2[abs(temp2$COR) > 0.25,  ]
#     num_sig_cors <- c(num_sig_cors, length(temp3[[1]]))
#     num_sig_cors_d[[k]] <- igraph::graph_from_data_frame(temp3[,c(1,2)], directed=FALSE)
#     background_as <- igraph::union(igraph::graph_from_data_frame(temp3[,c(1,2)], directed=FALSE), background_as)
# }

# pdata <- data.frame(cancer=all_cancer, count=num_sig_cors)

# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity",position="stack")+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of event-gene pairs\nthat significantly correlate", limits=c(0,max(pdata$count)+50)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')#guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#
# ggsave(p,filename=paste0("../results/Event_gene_pairs_cor.png"),width=3.5, height=3, dpi=400)


# ##--- figure for the overlap ---------------------------
# ##--- Event-gene pairs occuring in multiple cancer types --------------------------------
# combs <- list() ## store all combinations
# for(k in 1:length(all_cancer)){
#     combs[[k]] <- combn(all_cancer, k)
# }

# num_events_combs <- list()
# num_events_combs_counts <- c()

# for(k in 1:length(all_cancer)){

#     temp_comb <- as.data.frame(combs[[k]])
#     loop1 <- length(temp_comb)
#     loop2 <- length(temp_comb[[1]])
#     temp_unn <- igraph::make_empty_graph(directed=FALSE)
#     for(i in 1:loop1){
#         temp_ovl <- num_sig_cors_d[[which(all_cancer == temp_comb[[i]][1])]]
#         if(loop2 > 1){
#             for(j in 2:loop2){
#                 temp_ovl <- igraph::intersection(temp_ovl, num_sig_cors_d[[which(all_cancer == temp_comb[[i]][j])]])
#             }
#         }
#         temp_unn <- igraph::union(temp_unn, temp_ovl)
#     }

#     num_events_combs[[k]] <- temp_unn
#     num_events_combs_counts <- c(num_events_combs_counts, igraph::ecount(temp_unn))
#     cat('Cancer',k,'of',length(all_cancer),'done\n')
# }

# pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)

# p <- ggplot(pdata, aes(cancer, count, label=count)) + 
# geom_bar(stat="identity")+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# maxv <- max(pdata$count)
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="# of cancer types") + 
# scale_y_continuous(name="# of significantly correlated\nevent-gene pairs", limits=c(0,maxv+500)) +
# geom_text(position=position_dodge(0.9), hjust=0, vjust=0, angle=75, size=2.5)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
# ggsave(p,filename=paste0("../results/Event_gene_pairs_cor_overlap.png"),width=7, height=3, dpi=400)




