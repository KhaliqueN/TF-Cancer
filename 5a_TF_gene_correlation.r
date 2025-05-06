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
if(dir.exists(store_dir)){
    unlink(store_dir, recursive=TRUE)
}
dir.create(store_dir)

fdr <- 0.05
diff <- 0.1
rnd_expr <- 1000
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
tf_ensemb_map <- as.data.frame(data.table::fread('../data/TF_ensembl_uniprot.txt', sep='\t'))
uniprot_ensembl_map <- data.table::fread('../data/ensembl_name_map.txt')
symbol_ensembl_map <- unique(uniprot_ensembl_map[,c(1,8)])
symbol_ensembl_map <- symbol_ensembl_map[which(symbol_ensembl_map$HGNC_symbol != ''),]

all_filesxx <- gtools::mixedsort(list.files(psi_input, pattern='*filtered_PSI.txt', full.names=TRUE))
all_cancer <- substr(basename(all_filesxx), 1,4)

##-- contradictory event pairs ------------------------------
cont_pairs <- data.table::fread(paste0('../data/contradictory_event_pair_DBD.txt'))


all_files_genes <- gtools::mixedsort(list.files('../../../public_data/TCGA_normalized_counts', full.names=TRUE))
## Some of the event coordinates from TCGA splice seq does not overlap to any of the exons mapped from the canonical protein
## sequence from Uniprot. This is because the canonical seqeunce is not always the longest.

##--- grn without direction downloaded from here: https://tflink.net/download/
grn <- as.data.frame(data.table::fread('../data/TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv'))
grn_filt <- grn[grn$`Name.TF` %in% tfs$Gene_Symbol,]
##--- 14,439 edges ----
genes2consider <- union(grn_filt$`Name.TF`, grn_filt$`Name.Target`)

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_filesxx[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- as.data.frame(tempx[whx,])
    wh <- which(colnames(tempy) %like% '_Norm')
    tempy1 <- tempy[, -wh]
    wh <- which(colnames(tempy1) %like% 'TCGA')

    ## compute delta psi values with respect to the median normal psi value
    for(j in 1:length(wh)){
        txc <- tempy1[[wh[j]]]
        whx <- which(txc != 'null')
        txc[whx] <- as.numeric(txc[whx]) - tempy1$MEDIAN_NORMAL[whx]
        tempy1[wh[j]] <- txc
    }

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

    ##--- correlation between the delta psi values and the gene expressions -----
    tempy3x <- as.data.frame(sapply(tempy3,as.numeric))
    tcount <- apply(tempy3, 1, function(x) length(which(!is.na(as.numeric(x))))) ## count the number of non-NAs --> Should be at least 10 as
    ## I computed the differential PSI values for only those events that have at elast 10 non-NA values --> see script "1b_TCGA_data.r"

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

        # cat('Event',j,'of',length(tempy1[[1]]), 'done\n')
    }

    tempd <- data.frame(ASID=temp_as, GENE=temp_gene, COR=temp_cor, CPVAL=temp_pval, RPVAL=temp_rnd/rnd_expr)
    tempd <- tempd[complete.cases(tempd),]
    tempd$CFDR <- p.adjust(tempd$CPVAL, 'fdr')
    tempd$RFDR <- p.adjust(tempd$RPVAL, 'fdr')

    # wh1 <- which(tempd$RFDR < fdr)
    # wh2 <- which(tempd$CFDR < fdr)
    # wh <- intersect(wh1, wh2)
    # tempdf <- tempd[wh, ]
    # tempdf <- tempdf[order(-abs(tempdf$COR)),]

    data.table::fwrite(tempd, paste0(store_dir,'/',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)
    cat('Cancer',k, 'of',length(all_cancer),'done\n')

}


# ### *****  Number of event-gene pairs **********************************************************************
# # *****************************************************************************
# # *****************************************************************
# # ***************************************************

##--------------------
all_sigs <- gtools::mixedsort(list.files(store_dir, pattern='*.txt', full.names=TRUE))
num_sig_cors <- c()
num_sig_cors_d <- list()
background_as <- igraph::make_empty_graph(n = 0, directed = FALSE)

for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_sigs[k], sep='\t')
    temp1 <- temp[temp$RFDR < 0.05, ]
    temp2 <- temp1[temp1$CFDR < 0.05, ]
    num_sig_cors <- c(num_sig_cors, length(temp2[[1]]))
    num_sig_cors_d[[k]] <- igraph::graph_from_data_frame(temp2[,c(1,2)], directed=FALSE)
    background_as <- igraph::union(igraph::graph_from_data_frame(temp2[,c(1,2)], directed=FALSE), background_as)
}

pdata <- data.frame(cancer=all_cancer, count=num_sig_cors)

p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of event-gene pairs\nthat significantly correlate", limits=c(0,max(pdata$count)+250)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#
ggsave(p,filename=paste0("../results/Event_gene_pairs_cor.png"),width=3.5, height=3, dpi=400)


##--- figure for the overlap ---------------------------
##--- Event-gene pairs occuring in multiple cancer types --------------------------------
combs <- list() ## store all combinations
for(k in 1:length(all_cancer)){
    combs[[k]] <- combn(all_cancer, k)
}

num_events_combs <- list()
num_events_combs_counts <- c()

for(k in 1:length(all_cancer)){

    temp_comb <- as.data.frame(combs[[k]])
    loop1 <- length(temp_comb)
    loop2 <- length(temp_comb[[1]])
    temp_unn <- igraph::make_empty_graph(directed=FALSE)
    for(i in 1:loop1){
        temp_ovl <- num_sig_cors_d[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- igraph::intersection(temp_ovl, num_sig_cors_d[[which(all_cancer == temp_comb[[i]][j])]])
            }
        }
        temp_unn <- igraph::union(temp_unn, temp_ovl)
    }

    num_events_combs[[k]] <- temp_unn
    num_events_combs_counts <- c(num_events_combs_counts, igraph::ecount(temp_unn))
    cat('Cancer',k,'of',length(all_cancer),'done\n')
}

pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)

p <- ggplot(pdata, aes(cancer, count, label=count)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(pdata$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of splicing events", limits=c(0,maxv+500)) +
geom_text(position=position_dodge(0.9), hjust=0, vjust=0, angle=75, size=2.5)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0("../results/Event_gene_pairs_cor_overlap.png"),width=7, height=3, dpi=400)


### ******************************************************
### (2) Opposite effects on the gene ---
ASID1 <- c()
ASID2 <- c()
TF_SYMBOL <- c()
TARGET_SYMBOL <- c()
COR1 <- c()
COR2 <- c()
TCANCER <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_sigs[k])
    temp <- temp[temp$RFDR < 0.05, ]
    temp <- temp[temp$CFDR < 0.05, ]
    tempsig <- data.table::fread(all_filesxx[k], sep='\t')
    tempsig <- tempsig[tempsig$as_id %in% temp$ASID, ]
    temp_tfs <- c()
    for(j in 1:length(temp$ASID)){
        temp_tfs <- c(temp_tfs, tempsig[tempsig$as_id == temp$ASID[j], ]$symbol)
    }
    temp$TF <- temp_tfs
    tfg <- unique(temp$TF)

    for(j in 1:length(tfg)){

        temp1 <- temp[temp$TF == tfg[j], ]
        temp1a <- temp1[temp1$COR > 0, ]
        temp1b <- temp1[temp1$COR < 0, ]
        trg <- unique(temp1a$GENE)

        for(i in 1:length(trg)){
            t1 <- temp1a[temp1a$GENE == trg[i], ]
            t2 <- temp1b[temp1b$GENE == trg[i], ]
            if(nrow(t1) != 0 & nrow(t2) != 0){
                for(m in 1:length(t1[[1]])){
                    for(n in 1:length(t2[[1]])){
                        ASID1 <- c(ASID1, t1$ASID[m])
                        ASID2 <- c(ASID2, t2$ASID[n])
                        TF_SYMBOL <- c(TF_SYMBOL, t1$TF[m])
                        TARGET_SYMBOL <- c(TARGET_SYMBOL, t1$GENE[m])
                        TCANCER <- c(TCANCER, all_cancer[k])
                        COR1 <- c(COR1, t1$COR[m])
                        COR2 <- c(COR2, t2$COR[n])
                    }
                }
            }
        }
    }

    cat('Cancer',k, 'of',length(all_cancer),'done\n')

}

ydata_diff <- data.frame(CANCER=TCANCER, TF=TF_SYMBOL, TARGET=TARGET_SYMBOL, ASID1=ASID1, ASID2=ASID2, COR1=COR1, COR2=COR2)

data.table::fwrite(ydata_diff, paste0('../data/contradictory_event_pair-gene_comb.txt'), row.names=FALSE, quote=FALSE, sep='\t')

##-- correlation plot ----
##--- scatter plot ----
ydata <- ydata_diff
p <- ggplot(ydata, aes(COR2, COR1, color=CANCER)) + 
geom_point(alpha=0.6)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="TF PSI correlation with \ngene expression patterns", limits=c(-1,0)) +
scale_y_continuous(name="TF PSI correlation with \ngene expression patterns", limits=c(0,1)) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Event_gene_pairs_cor_cors.png"),width=5, height=3, dpi=400)


##-- number of contradictory event pairs ---
pdata <- plyr::count(ydata_diff[[1]])
p <- ggplot(pdata, aes(x, freq)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of contradictory event pairs", limits=c(0,max(pdata$freq)+50)) +
geom_text(aes(label=freq), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0("../results/Event_gene_pairs_cor_cont.png"),width=3.5, height=3, dpi=400)


##--- Pairs of event-gene pairs affecting the same gene that are contradictory in terms of 
##--- splicing ----
fdatay <- ydata_diff
fdatay$X1 <- paste0(fdatay$TF, '_', fdatay$ASID1)
fdatay$X2 <- paste0(fdatay$TF, '_', fdatay$ASID2)
cont_pairs$X1 <- unlist(lapply(cont_pairs$AS1, function(x) substr(x, 1, nchar(x)-3)))
cont_pairs$X2 <- unlist(lapply(cont_pairs$AS2, function(x) substr(x, 1, nchar(x)-3)))

num_sig_cors_d <- list()
# num_sig_cors_d_count <- c()
for(k in 1:length(all_cancer)){
    fdatayt <- fdatay[fdatay$CANCER == all_cancer[k], ]
    fdatayu <- cont_pairs[cont_pairs$CANCER == all_cancer[k], ]
    fdatayt_gr <- igraph::graph_from_data_frame(fdatayt[,c(8,9)], directed=FALSE)
    fdatayu_gr <- igraph::graph_from_data_frame(fdatayu[,c(8,9)], directed=FALSE)
    num_sig_cors_d[[k]] <- igraph::intersection(fdatayt_gr, fdatayu_gr)
    # num_sig_cors_d_count <- c(num_sig_cors_d_count, igraph::ecount(igraph::intersection(fdatayt_gr, fdatayu_gr)))
}


num_sig_cors_dx <- list()
num_sig_cors_dx_c <- c()
for(k in 1:length(all_cancer)){
    temp <- igraph::as_data_frame(num_sig_cors_d[[k]])
    temp1 <- fdatay[fdatay$CANCER == all_cancer[k], ]
    pos <- c()
    for(j in 1:length(temp[[1]])){
        wh1 <- which(temp1$X1 == temp[[1]][j])
        wh2 <- which(temp1$X2 == temp[[2]][j])
        wha <- intersect(wh1, wh2)
        wh1 <- which(temp1$X2 == temp[[1]][j])
        wh2 <- which(temp1$X1 == temp[[2]][j])
        whb <- intersect(wh1, wh2)
        wh <- union(wha, whb)
        pos <- c(pos, wh)
    }

    temp_gr <- temp1[pos, ]
    temp_gr$e1 <- paste0(temp_gr$ASID1,'_',temp_gr$TARGET)
    temp_gr$e2 <- paste0(temp_gr$ASID2,'_',temp_gr$TARGET)
    num_sig_cors_dx[[k]] <- igraph::graph_from_data_frame(temp_gr[,c(10,11)], directed=FALSE)
    num_sig_cors_dx_c <- c(num_sig_cors_dx_c, length(temp_gr[[1]]))
}


##-- number of contradictory event pairs ---
pdata <- data.frame(x=all_cancer, freq=num_sig_cors_dx_c)#plyr::count(ydata_diff[[1]])
p <- ggplot(pdata, aes(x, freq)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of contradictory event pairs showing \ncontradictory event-gene associations", limits=c(0,max(pdata$freq)+50)) +
geom_text(aes(label=freq), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0("../results/Contradictory_event_gene_pairs_cor_cont.png"),width=3.5, height=3, dpi=400)


##---- overlap ----
##--- Splicing events occuring in multiple cancer types --------------------------------
# pdatax <- ydata_diff
# pdatax$e1 <- paste0(pdatax$ASID1,'_',pdatax$TARGET)
# pdatax$e2 <- paste0(pdatax$ASID2,'_',pdatax$TARGET)

combs <- list() ## store all combinations
for(k in 1:length(all_cancer)){
    combs[[k]] <- combn(all_cancer, k)
}

# ##-- create undirected graphs --
# library(igraph)
# cancer_gr <- list()
# for(k in 1:length(all_cancer)){
#     temp <- pdatax[pdatax$CANCER == all_cancer[k],]
#     temp1 <- temp[,c(8,9)]
#     cancer_gr[[k]] <- igraph::graph_from_data_frame(temp1, directed=FALSE)
# }

cancer_gr <- num_sig_cors_dx
num_events_combs <- list()
num_events_combs_counts <- c()
for(k in 1:length(all_cancer)){

    temp_comb <- as.data.frame(combs[[k]])
    loop1 <- length(temp_comb)
    loop2 <- length(temp_comb[[1]])
    temp_unn <- igraph::make_empty_graph(n = 0, directed = FALSE)

    for(i in 1:loop1){
        temp_ovl <- cancer_gr[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- igraph::intersection(temp_ovl, cancer_gr[[which(all_cancer == temp_comb[[i]][j])]], keep.all.vertices = TRUE)
            }
        }
        temp_unn <- igraph::union(temp_unn, temp_ovl)
    }

    num_events_combs[[k]] <- temp_unn
    num_events_combs_counts <- c(num_events_combs_counts, igraph::ecount(temp_unn))
    cat('Cancer',k,'of',length(all_cancer),'done\n')
}

pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of contradictory event pairs", limits=c(0,max(pdata$count)+220)) +
geom_text(aes(label=count), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0("../results/Contradictory_event_gene_pairs_cor_cont_overlap.png"),width=7, height=3, dpi=400)


