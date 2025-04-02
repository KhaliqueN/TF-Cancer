##############################################################################################
# Purpose: Correlation between patient-specific delta PSI values of TFs and patient-specific 
# gene expression patterns
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)
library(dorothea)
library(pheatmap)

save_dir <- '../results'
psi_input <- '../data/PSI_data'

## PCA function -----------------
topca <- function(inputfile, variationThreshold){

  x <- inputfile
  noRows <- nrow(x)
  noColumns <- ncol(x)
  x.pca <- prcomp(x, scale.=TRUE)
  sumEigen <- sum(x.pca$sdev^2)
  sumEigenCurrent <- 0
  for(i in 1:noColumns)
  {
    sumEigenCurrent <- sumEigenCurrent + x.pca$sdev[i]^2
    if(sumEigenCurrent/sumEigen >= variationThreshold)
    {
      break
    }
  }
  noColumnsReduced <- i
  if(noColumnsReduced < 2) {
    noColumnsReduced <- 2
    # cat("No of reduced dimensions is corrected to 2 to ensure cosine similarity-based comparison\n")
  }
  
  return(x.pca$x[,1:noColumnsReduced])

}


fdr <- 0.05
diff <- 0.1
rnd_expr <- 1000
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
tf_ensemb_map <- as.data.frame(data.table::fread('../data/TF_ensembl_uniprot.txt', sep='\t'))

all_filesxx <- gtools::mixedsort(list.files(psi_input, pattern='*filtered_PSI.txt', full.names=TRUE))
all_cancer <- substr(basename(all_filesxx), 1,4)

all_files_genes <- gtools::mixedsort(list.files('../../../public_data/TCGA_normalized_counts', full.names=TRUE))
## Some of the event coordinates from TCGA splice seq does not overlap to any of the exons mapped from the canonical protein
## sequence from Uniprot. This is because the canonical seqeunce is not always the longest.

CORRS <- c()
# Z_CORRS <- c()
# P_VALS <- c()
CORRS_VAL <- c()
EVENTS <- c()
CANCER <- c()
RND_SIG <- c()

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

    ## gene expression ----
    temp <- as.data.frame(data.table::fread(all_files_genes[k], sep='\t'))
    labs <- temp[length(temp)]
    temp1 <- cbind(labs,temp[,-length(temp)])
    tgenes <- unlist(lapply(strsplit(temp1[[1]],'[.]'),'[[',1))
    temp <- temp1[,-1]
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

    ## compute PCA of the gene expression matrix --------------
    cancer_temp_t <- t(cancer_temp_filt)
    ## remove zero deviation columns
    allsds <- matrixStats::colSds(cancer_temp_t)
    whm <- which(allsds != 0)
    tgenes <- tgenes[whm]
    cancer_temp_t <- cancer_temp_t[,whm]
    gene_expr_pca <- as.data.frame(topca(cancer_temp_t, 0.9))


    ##--- correlation between the delta psi values and the first PCA component
    all_cors <- c()
    all_cors_vals <- c()
    # z_cors <- c()
    # p_vals <- c()
    rnd_cors_est <- c()

    for(j in 1:length(tempy3[[1]])){
        txc <- tempy3[j,]
        wh <- which(txc != 'null')
        txc <- as.numeric(txc[wh])
        temp_cor <- cor.test(x=gene_expr_pca[[1]][wh], y=txc, method = 'spearman')
        
        if(is.na(temp_cor$p.value)){
            all_cors <- c(all_cors, 1)
            all_cors_vals <- c(all_cors_vals, 0)
        }else{
            all_cors <- c(all_cors, temp_cor$p.value)
            all_cors_vals <- c(all_cors_vals, temp_cor$estimate)
        }

        ##-- compute cor on randomized data
        # t_rnd_cors <- c()
        t_rnd_cors_c <- c()

        for(i in 1:rnd_expr){
            temp_cor_rnd <- cor.test(x=sample(gene_expr_pca[[1]][wh]), y=sample(txc), method = 'spearman')
            # t_rnd_cors <- c(t_rnd_cors, temp_cor_rnd$estimate)
            t_rnd_cors_c <- c(t_rnd_cors_c, temp_cor_rnd$p.value)
        }

        # rnd_sd <- sd(t_rnd_cors)
        # rnd_mn <- mean(t_rnd_cors)
        # zScore <- (temp_cor$estimate-rnd_mn)/rnd_sd
        # z_cors <- c(z_cors, zScore)
        # p_vals <- c(p_vals, 2*pnorm(abs(zScore), mean = rnd_mn, sd = rnd_sd, lower.tail = FALSE))
        rnd_cors_est <- c(rnd_cors_est, length(which(t_rnd_cors_c <= temp_cor$p.value)))

        cat('Cancer',k,':',all_cancer[k],': Gene',j, 'of',length(tempy3[[1]]),'done\n')
    }

    CORRS <- c(CORRS, all_cors)
    # Z_CORRS <- c(Z_CORRS, z_cors)
    # P_VALS <- c(P_VALS, p_vals)
    CORRS_VAL <- c(CORRS_VAL, all_cors_vals)
    EVENTS <- c(EVENTS, tempy1$as_id)
    CANCER <- c(CANCER, rep(all_cancer[k], length(all_cors)))
    RND_SIG <- c(RND_SIG, rnd_cors_est)

}

# pdata <- data.frame(cancer=CANCER, Event=EVENTS, Zscore=Z_CORRS, Zpvalue=P_VALS, Spearman_value=CORRS_VAL, Spearman_pvalue=CORRS, Random_sig=RND_SIG/rnd_expr)
pdata <- data.frame(cancer=CANCER, Event=EVENTS, Spearman_value=CORRS_VAL, Spearman_pvalue=CORRS, Random_sig=RND_SIG/rnd_expr)

data.table::fwrite(pdata, '../data/TF_deltaPSI_geneExpr_cor_1000.txt', sep='\t', row.names=FALSE, quote=FALSE)

pdata <- data.table::fread('../data/TF_deltaPSI_geneExpr_cor_1000.txt')
## FDR correction -------
pdata$SPEARMAN_FDR <- p.adjust(pdata$Spearman_pvalue, 'fdr')
# pdata$Z_FDR <- p.adjust(pdata$Zpvalue, 'fdr')
pdata$RND_FDR <- p.adjust(pdata$Random_sig, 'fdr')

pdatax <- pdata[pdata$SPEARMAN_FDR < 0.05, ]
pdatay <- pdatax[pdatax$RND_FDR < 0.05, ]



##---- TF splicing events that show high correlation with gene expression patterns ----
num_sig_events_pos <- c()
num_sig_events_neg <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_filesxx[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]

    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy <- tempy[tempy$as_id %in% pdatay[pdatay$cancer == all_cancer[k], ]$Event, ]

    tempy_pos <- tempy[tempy$MEDIAN_DIFF > 0, ]
    tempy_neg <- tempy[tempy$MEDIAN_DIFF < 0, ]

    num_sig_events_pos <- c(num_sig_events_pos, length(tempy_pos[[1]]))
    num_sig_events_neg <- c(num_sig_events_neg, length(tempy_neg[[1]]))

}

##--- plot the number of TF splicing events significantly correlated with gene expression patterns ------------
pdata1 <- data.frame(cancer=all_cancer, count=num_sig_events_pos, flag=rep('High in cancer'))
pdata2 <- data.frame(cancer=all_cancer, count=num_sig_events_neg, flag=rep('Low in cancer'))
pdata <- rbind(pdata1,pdata2)

p <- ggplot(pdata, aes(cancer, count, fill=flag, label=count)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(pdata1$count+pdata2$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of significant splicing events", limits=c(0,maxv)) +
geom_text(position=position_stack(vjust=0.5), size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Sig_events_TF_cor.png"),width=5.5, height=3, dpi=400)


##-------- Pairwise overlap of slicing events affecting TFs across cancer types --------------------------------------------
TF_splicing_events_pos <- list()
TF_splicing_events_neg <- list()
background_as_pos <- c()
background_as_neg <- c() 

for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_filesxx[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]

    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy <- tempy[tempy$as_id %in% pdatay[pdatay$cancer == all_cancer[k], ]$Event, ]

    tempy_pos <- tempy[tempy$MEDIAN_DIFF > 0,]
    tempy_neg <- tempy[tempy$MEDIAN_DIFF < 0,]

    TF_splicing_events_pos[[k]] <- tempy_pos$as_id
    background_as_pos <- union(background_as_pos, tempy_pos$as_id)

    TF_splicing_events_neg[[k]] <- tempy_neg$as_id
    background_as_neg <- union(background_as_neg, tempy_neg$as_id)

}


##--- figure for the overlap ---------------------------
c1 <- c()
c2 <- c()
oval <- c()
pval <- c()
pvalu <- c()
loop1 <- length(all_cancer)-1
loop2 <- length(all_cancer)

for(ov1 in 1:loop1){
    m <- ov1+1
    for(ov2 in m:loop2){
        e1 <- length(TF_splicing_events_pos[[ov1]])
        e2 <- length(TF_splicing_events_pos[[ov2]])
        c1 <- c(c1, paste0(all_cancer[ov1],' (', e1,')'))
        c2 <- c(c2, paste0(all_cancer[ov2],' (', e2,')'))
        inter_e <- intersect(TF_splicing_events_pos[[ov1]], TF_splicing_events_pos[[ov2]])
        oval <- c(oval, length(inter_e))

        hyp <- phyper(length(inter_e)-1,e1,length(background_as_pos)-e1,e2,lower.tail = FALSE)
        hyp <- signif(hyp, digits=3)
        pval <- c(pval, hyp)

        ## under representation
        hyp <- phyper(length(inter_e),e1,length(background_as_pos)-e1,e2,lower.tail = TRUE)
        hyp <- signif(hyp, digits=3)
        pvalu <- c(pvalu, hyp)
    # }
    }
}

qval <- p.adjust(pval, 'fdr')
qvalu <- p.adjust(pvalu, 'fdr')
qvalx <- rep('Expected by chance',length(pval))
for(i in 1:length(pval)){
    if(qval[i] < fdr){
        qvalx[i] <- 'Significantly \nmore overlapping'
    }
    if(qvalu[i] < fdr){
        qvalx[i] <- 'Significantly \nless overlapping'
    }
}


pdata <- data.frame(c1=c1, c2=c2, val=oval, pval=qvalx, qval=qval)

cols <- c('#377eb8','#4daf4a','#e41a1c')#rev(brewer.pal(3,"Spectral"))
p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pval),colour = "white")+
  theme(legend.text=element_text(size=8))+scale_fill_manual(values=cols,drop=FALSE)
basesize <- 8
p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type", y = "Cancer type") +
  scale_y_discrete() +
  scale_x_discrete()+coord_fixed()+
  guides(fill=guide_legend(title=""))+
  geom_text(data=pdata,aes(y=c2,label=val),size=1.5)+
  theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0.5, colour = "black"))
  # +guides(fill='none')
ggsave(p,filename=paste0(save_dir,'/Sig_events_TFs_overlap_pos_cor.png'),width=4, height=2.5, dpi=400)


c1 <- c()
c2 <- c()
oval <- c()
pval <- c()
pvalu <- c()
loop1 <- length(all_cancer)-1
loop2 <- length(all_cancer)

for(ov1 in 1:loop1){
    m <- ov1+1
    for(ov2 in m:loop2){
        e1 <- length(TF_splicing_events_neg[[ov1]])
        e2 <- length(TF_splicing_events_neg[[ov2]])
        c1 <- c(c1, paste0(all_cancer[ov1],' (', e1,')'))
        c2 <- c(c2, paste0(all_cancer[ov2],' (', e2,')'))
        inter_e <- intersect(TF_splicing_events_neg[[ov1]], TF_splicing_events_neg[[ov2]])
        oval <- c(oval, length(inter_e))

        hyp <- phyper(length(inter_e)-1,e1,length(background_as_neg)-e1,e2,lower.tail = FALSE)
        hyp <- signif(hyp, digits=3)
        pval <- c(pval, hyp)

        ## under representation
        hyp <- phyper(length(inter_e),e1,length(background_as_neg)-e1,e2,lower.tail = TRUE)
        hyp <- signif(hyp, digits=3)
        pvalu <- c(pvalu, hyp)
    # }
    }
}

qval <- p.adjust(pval, 'fdr')
qvalu <- p.adjust(pvalu, 'fdr')
qvalx <- rep('Expected by chance',length(pval))
for(i in 1:length(pval)){
    if(qval[i] < fdr){
        qvalx[i] <- 'Significantly \nmore overlapping'
    }
    if(qvalu[i] < fdr){
        qvalx[i] <- 'Significantly \nless overlapping'
    }
}


pdata <- data.frame(c1=c1, c2=c2, val=oval, pval=qvalx, qval=qval)

cols <- c('#377eb8','#e41a1c')#rev(brewer.pal(3,"Spectral"))
p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pval),colour = "white")+
  theme(legend.text=element_text(size=8))+scale_fill_manual(values=cols,drop=FALSE)
basesize <- 8
p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type", y = "Cancer type") +
  scale_y_discrete() +
  scale_x_discrete()+coord_fixed()+
  guides(fill=guide_legend(title=""))+
  geom_text(data=pdata,aes(y=c2,label=val),size=1.5)+
  theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0.5, colour = "black"))
  # +guides(fill='none')
ggsave(p,filename=paste0(save_dir,'/Sig_events_TFs_overlap_neg_cor.png'),width=4, height=2.5, dpi=400)


c1 <- c()
c2 <- c()
oval <- c()
pval <- c()
pvalu <- c()
loop1 <- length(all_cancer)-1
loop2 <- length(all_cancer)
background_as <- union(background_as_neg, background_as_pos)
for(ov1 in 1:loop1){
    m <- ov1+1
    for(ov2 in m:loop2){
        e1 <- length(TF_splicing_events_pos[[ov1]])
        e2 <- length(TF_splicing_events_neg[[ov2]])
        c1 <- c(c1, paste0(all_cancer[ov1],' (', e1,')'))
        c2 <- c(c2, paste0(all_cancer[ov2],' (', e2,')'))
        inter_e <- intersect(TF_splicing_events_pos[[ov1]], TF_splicing_events_neg[[ov2]])
        oval <- c(oval, length(inter_e))

        hyp <- phyper(length(inter_e)-1,e1,length(background_as)-e1,e2,lower.tail = FALSE)
        hyp <- signif(hyp, digits=3)
        pval <- c(pval, hyp)

        ## under representation
        hyp <- phyper(length(inter_e),e1,length(background_as)-e1,e2,lower.tail = TRUE)
        hyp <- signif(hyp, digits=3)
        pvalu <- c(pvalu, hyp)
    # }
    }
}

qval <- p.adjust(pval, 'fdr')
qvalu <- p.adjust(pvalu, 'fdr')
qvalx <- rep('Expected by chance',length(pval))
for(i in 1:length(pval)){
    if(qval[i] < fdr){
        qvalx[i] <- 'Significantly \nmore overlapping'
    }
    if(qvalu[i] < fdr){
        qvalx[i] <- 'Significantly \nless overlapping'
    }
}


pdata <- data.frame(c1=c1, c2=c2, val=oval, pval=qvalx, qval=qval)

cols <- c('#377eb8','#4daf4a','#e41a1c')#rev(brewer.pal(3,"Spectral"))
p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pval),colour = "white")+
  theme(legend.text=element_text(size=8))+scale_fill_manual(values=cols,drop=FALSE)
basesize <- 8
p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type (high PSI)", y = "Cancer type (low PSI)") +
  scale_y_discrete() +
  scale_x_discrete()+coord_fixed()+
  guides(fill=guide_legend(title=""))+
  geom_text(data=pdata,aes(y=c2,label=val),size=1.5)+
  theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0.5, colour = "black"))
  # +guides(fill='none')
ggsave(p,filename=paste0(save_dir,'/Sig_events_TFs_overlap_pos_neg_cor.png'),width=4, height=2.5, dpi=400)



##--- Splicing events occuring in multiple cancer types --------------------------------
combs <- list() ## store all combinations
for(k in 1:length(all_cancer)){
    combs[[k]] <- combn(all_cancer, k)
}

##-- speed up --read all files in memory!!
all_files_store <- list()
for(k in 1:length(all_files)){
    temp <- as.data.frame(data.table::fread(all_filesxx[k], sep='\t'))
    ll <- length(temp) - 3
    all_files_store[[k]] <- read.table(all_filesxx[k], colClasses = c('NULL','character','character',rep('NULL', ll)), header = TRUE, fill=TRUE)
}

num_events_combs_pos <- list()
num_events_combs_counts <- c()

for(k in 1:length(all_cancer)){

    temp_comb <- as.data.frame(combs[[k]])
    loop1 <- length(temp_comb)
    loop2 <- length(temp_comb[[1]])
    temp_unn <- c()
    for(i in 1:loop1){
        temp_ovl <- TF_splicing_events_pos[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- intersect(temp_ovl, TF_splicing_events_pos[[which(all_cancer == temp_comb[[i]][j])]])
            }
        }
        temp_unn <- union(temp_unn, temp_ovl)
    }
    num_events_combs_pos[[k]] <- temp_unn
    num_events_combs_counts <- c(num_events_combs_counts, length(temp_unn))
}

pdata1 <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)

num_events_combs_neg <- list()
num_events_combs_counts <- c()

for(k in 1:length(all_cancer)){

    temp_comb <- as.data.frame(combs[[k]])
    loop1 <- length(temp_comb)
    loop2 <- length(temp_comb[[1]])
    temp_unn <- c()
    for(i in 1:loop1){
        temp_ovl <- TF_splicing_events_neg[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- intersect(temp_ovl, TF_splicing_events_neg[[which(all_cancer == temp_comb[[i]][j])]])
            }
        }
        temp_unn <- union(temp_unn, temp_ovl)
    }
    num_events_combs_neg[[k]] <- temp_unn
    num_events_combs_counts <- c(num_events_combs_counts, length(temp_unn))
}

pdata2 <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)
pdata1$flag <- rep('High in cancer', length(pdata1[[1]]))
pdata2$flag <- rep('Low in cancer', length(pdata2[[1]]))
pdata <- rbind(pdata1,pdata2)
p <- ggplot(pdata, aes(cancer, count, fill=flag, label=count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(max(pdata1$count), max(pdata2$count))
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of splicing events", limits=c(0,maxv+100)) +
geom_text(position=position_dodge(0.9), hjust=0, vjust=0, angle=75, size=2.5)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Events_perturbing_overlap_cor.png"),width=7, height=3, dpi=400)


##----- PSI profile of the pancancer (at least 10 cancer types) events -------------
events_to_consider <- num_events_combs_pos[[10]]
tcancer <- c()
tvalue <- c()
tevent <- c()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_filesxx[k], sep='\t')
    temp1 <- temp[temp$as_id %in% events_to_consider, ]
    temp1$ID <- paste0(temp1$symbol,'_',temp1$as_id,'_',temp1$splice_type)
    tcancer <- c(tcancer, rep(all_cancer[k], length(temp1[[1]])))
    tvalue <- c(tvalue, temp1$MEDIAN_DIFF)
    tevent <- c(tevent, temp1$ID)
}

pdata1 <- data.frame(tcancer=tcancer, val=tvalue, event=tevent)

events_to_consider <- num_events_combs_neg[[10]]
tcancer <- c()
tvalue <- c()
tevent <- c()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_filesxx[k], sep='\t')
    temp1 <- temp[temp$as_id %in% events_to_consider, ]
    temp1$ID <- paste0(temp1$symbol,'_',temp1$as_id,'_',temp1$splice_type)
    tcancer <- c(tcancer, rep(all_cancer[k], length(temp1[[1]])))
    tvalue <- c(tvalue, temp1$MEDIAN_DIFF)
    tevent <- c(tevent, temp1$ID)
}

pdata2 <- data.frame(tcancer=tcancer, val=tvalue, event=tevent)

pdata1$flag <- rep('Positive', length(pdata1[[1]]))
pdata2$flag <- rep('Negative', length(pdata2[[1]]))
pdata <- rbind(pdata1, pdata2)

pdata_mat <- as.data.frame(matrix(nrow=length(unique(pdata$event)), ncol=length(all_cancer),0))
rownames(pdata_mat) <- gtools::mixedsort(unique(pdata$event))
colnames(pdata_mat) <- gtools::mixedsort(unique(all_cancer))

for(k in 1:length(pdata[[1]])){
    wh1 <- which(rownames(pdata_mat) == pdata$event[k])
    wh2 <- which(colnames(pdata_mat) == pdata$tcancer[k])
    pdata_mat[wh1, wh2] <- as.numeric(pdata$val[k])
}

# pdata_mat <- cbind(data.frame(Event=unique(pdata$event)), pdata_mat)
pdx <- as.matrix(pdata_mat)

p <- pheatmap::pheatmap(pdx)
ggsave(p,filename=paste0(save_dir, "/Pancancer_events_cor.png"),width=7, height=5, dpi=400)


##----- unique splicing events -----------------------------------------------------
num_events_unq <- list()
num_events_unq_counts <- c()
unique_events_tf <- list()

for(k in 1:length(all_cancer)){
    temp_unq <- TF_splicing_events_pos[[k]]
    temp_un <- c()
    for(j in 1:length(all_cancer)){
        if(j != k){
            temp_un <- union(temp_un, TF_splicing_events_pos[[j]])
        }
    }
    temp_unq <- setdiff(temp_unq, temp_un)
    num_events_unq[[k]] <- temp_unq
    num_events_unq_counts <- c(num_events_unq_counts, length(temp_unq))

}

pdata1 <- data.frame(cancer=all_cancer, count=num_events_unq_counts)

num_events_unq <- list()
num_events_unq_counts <- c()
unique_events_tf <- list()

for(k in 1:length(all_cancer)){
    temp_unq <- TF_splicing_events_neg[[k]]
    temp_un <- c()
    for(j in 1:length(all_cancer)){
        if(j != k){
            temp_un <- union(temp_un, TF_splicing_events_neg[[j]])
        }
    }
    temp_unq <- setdiff(temp_unq, temp_un)
    num_events_unq[[k]] <- temp_unq
    num_events_unq_counts <- c(num_events_unq_counts, length(temp_unq))

}

pdata2 <- data.frame(cancer=all_cancer, count=num_events_unq_counts)
pdata1$flag <- rep('High in cancer', length(pdata1[[1]]))
pdata2$flag <- rep('Low in cancer', length(pdata2[[1]]))
pdata <- rbind(pdata1,pdata2)
p <- ggplot(pdata, aes(cancer, count, fill=flag, label=count)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(pdata1$count+pdata2$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of significant splicing events", limits=c(0,maxv)) +
geom_text(position=position_stack(vjust=0.5), size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Unique_sig_events_TFs_cor.png"),width=5.5, height=3, dpi=400)

###------------------------------------------------------------------------------------------------




###--- Pairs of events affecting the same TF that show opposite coorrelation with gene expression --
p1 <- c()
p2 <- c()
c1 <- c()
c2 <- c()
tcancer <- c()

for(k in 1:length(all_cancer)){

    tempx <- pdatay[pdatay$cancer == all_cancer[k], ]

    temp <- data.table::fread(all_filesxx[k], sep='\t')
    temp1 <- temp[temp$as_id %in% tempx$Event, ]

    all_sym <- unique(temp1$symbol)

    for(j in 1:length(all_sym)){

        tt1 <- temp1[temp1$symbol == all_sym[j], ]$as_id
        tempy <- tempx[tempx$Event %in% tt1, ]

        if(nrow(tempy) > 1){
            loop1 <- length(tempy[[1]])-1
            loop2 <- length(tempy[[1]])

            for(i in 1:loop1){
                tx1 <- temp1[temp1$as_id == tempy$Event[i], ]
                mm <- i+1
                for(ii in mm:loop2){
                    tx2 <- temp1[temp1$as_id == tempy$Event[ii], ]
                    p1 <- c(p1, paste0(tx1$symbol,'_',tx1$as_id,'_',tx1$splice_type))
                    p2 <- c(p2, paste0(tx2$symbol,'_',tx2$as_id,'_',tx2$splice_type))
                    c1 <- c(c1, tempy$Spearman_value[i])
                    c2 <- c(c2, tempy$Spearman_value[ii])
                    tcancer <- c(tcancer, all_cancer[k])
                }
            }
        }
    }
    cat('Cancer',k,'of',length(all_cancer),'done\n')
}


tempdata <- data.frame(CANCER=tcancer, AS1=p1, AS2=p2, COR1=c1, COR2=c2)
ttdata <- tempdata
##--- scatter plot ----
p <- ggplot(tempdata, aes(c1, c2, color=CANCER)) + 
geom_point()+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="TF PSI correlation with \ngene expression patterns", limits=c(-1,1)) +
scale_y_continuous(name="TF PSI correlation with \ngene expression patterns", limits=c(-1,1)) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0(save_dir,"/Event_pairs_cor_with_gene_expr.png"),width=7, height=5, dpi=400)



##--- which of the event pairs of the same gene that show opposite correlation with gene expressions also have affected DBDs-----
store_val <- function(y1){
    if(nrow(y1) == 0){
        flag <- 0
    }else{
        y2 <- y1[y1$DBD != '-', ]
        if(nrow(y2) == 0){
            flag <- 0
        }else{
            flag <- unique(y2$DBD)
        }
    }
    # if(length(flag) > 1){
    #     print(paste(flag,collapse=','))
    # }
    return(paste(flag,collapse=','))
}

odbd1 <- c()
odbd2 <- c()
tcancer <- c()
as1 <- c()
as2 <- c()
cor1 <- c()
cor2 <- c()
counter <- 0
counterx <- 0
## 52 entries missed because of no uniprot ensembl match
for(k in 1:length(all_cancer)){

    tempxx <- ttdata[ttdata$CANCER == all_cancer[k],]
    wh1 <- which(tempxx$COR1 < 0)
    wh2 <- which(tempxx$COR2 > 0)
    wha <- intersect(wh1, wh2)
    wh1 <- which(tempxx$COR1 > 0)
    wh2 <- which(tempxx$COR2 < 0)
    whb <- intersect(wh1, wh2)
    wh <- union(wha, whb)
    tempyy <- tempxx[wh, ]
    counterx <- counterx+length(tempyy[[1]])

    for(j in 1:length(tempyy[[1]])){
        txx <- strsplit(tempyy[j,]$AS1,'[_]')[[1]][1]
        teves <- c(strsplit(tempyy[j,]$AS1,'[_]')[[1]][2], strsplit(tempyy[j,]$AS2,'[_]')[[1]][2])
        tst <- c(strsplit(tempyy[j,]$AS1,'[_]')[[1]][3], strsplit(tempyy[j,]$AS2,'[_]')[[1]][3])
        tuni <- tf_ensemb_map[tf_ensemb_map$HGNC_symbol == txx, ]$Uniprotswissprot[1]
        if(is.na(tuni)){
            counter <- counter+1
            next
        } ## some gene symbols of TFs were not mapped to Ensembl

        # if("RI" %in% tst){
        #     next
        # }
        tmap <- as.data.frame(data.table::fread(paste0('../data/uniprot_Ensembl_Exon_map_DBD_AS/',tuni,'_',all_cancer[k],'.txt')))
        
        for(i in 1:length(teves)){
            if(tst[i] == 'AP'){
                y1 <- tmap[tmap$AP == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }   
            }else if(tst[i] == 'AT'){
                y1 <- tmap[tmap$AT == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }  
            }else if(tst[i] == 'ES'){
                y1 <- tmap[tmap$ES == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }  
            }else if(tst[i] == 'AD'){
                y1 <- tmap[tmap$AD == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }  
            }else if(tst[i] == 'AA'){
                y1 <- tmap[tmap$AA == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }  
            }else if(tst[i] == 'ME'){
                y1 <- tmap[tmap$ME == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }  
            }else{
                if(i == 1){
                    odbd1 <- c(odbd1, 0)
                }else{
                    odbd2 <- c(odbd2, 0)
                } 
            }
        }

        as1 <- c(as1, tempyy[j,]$AS1)
        as2 <- c(as2, tempyy[j,]$AS2)
        tcancer <- c(tcancer, all_cancer[k])
        cor1 <- c(cor1, tempyy[j,]$COR1)
        cor2 <- c(cor2, tempyy[j,]$COR2)
    }

}

tdata <- data.frame(CANCER=tcancer, AS1=as1, AS2=as2, DBD1=odbd1, DBD2=odbd2, COR1=cor1, COR2=cor2)
wh1 <- which(tdata$DBD1 != 0)
wh2 <- which(tdata$DBD2 != 0)
wh <- union(wh1, wh2)
tdatax <- tdata[wh,]
## 144 cases where at least one one event affects the DBD
## 140 of the 144 are those cases where one event affects the DBD and the other does not
wh1a <- setdiff(wh1, wh2)
wh2a <- setdiff(wh2, wh1)
wh <- union(wh1a, wh2a)
tdatay <- tdata[wh, ]

pdata <- plyr::count(tdatay[[1]])

p <- ggplot(pdata, aes(x, freq)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of contradictory event pairs", limits=c(0,20)) +
geom_text(aes(label=freq), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Contradictory_event_pairs.png"),width=3.5, height=3, dpi=400)


##--- Splicing events occuring in multiple cancer types --------------------------------
combs <- list() ## store all combinations
for(k in 1:length(all_cancer)){
    combs[[k]] <- combn(all_cancer, k)
}

##-- create undirected graphs --
library(igraph)
cancer_gr <- list()
for(k in 1:length(all_cancer)){
    temp <- tdatay[tdatay$CANCER == all_cancer[k],]
    temp1 <- temp[,c(2,3)]
    cancer_gr[[k]] <- igraph::graph_from_data_frame(temp1, directed=FALSE)
}

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
scale_y_continuous(name="# of contradictory event pairs", limits=c(0,150)) +
geom_text(aes(label=count), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Contradictory_event_pairs_overlap.png"),width=7, height=3, dpi=400)





# ##--- Do the events affecting the same gene occur simultaneosly?
# p1 <- c()
# p2 <- c()
# c1 <- c()
# c2 <- c()
# cor1 <- c()
# cor2 <- c()
# tcancer <- c()
# ttdata <- tempdata

# for(k in 1:length(all_cancer)){

#     temp <- ttdata[ttdata$CANCER == all_cancer[k],]
#     as1 <- unlist(lapply(strsplit(temp$AS1,'[_]'),'[[',2))
#     as2 <- unlist(lapply(strsplit(temp$AS2,'[_]'),'[[',2))

#     tempx <- as.data.frame(data.table::fread(all_filesxx[k], sep='\t'))
#     wh <- which(colnames(tempx) %like% '_Norm')
#     tempy <- tempx[, -wh]
#     tempy <- tempy[14:length(tempy[[1]]), ]
#     wh <- which(colnames(tempy) %like% 'TCGA')

#     for(j in 1:length(temp[[1]])){
#         t1 <- tempy[tempy$as_id == as1[j], wh][1,]
#         t2 <- tempy[tempy$as_id == as2[j], wh][1,]
#         wh1 <- which(t1 == 'null')
#         wh2 <- which(t2 == 'null')
#         whxx <- union(wh1, wh2)
#         t1 <- as.numeric(unname(unlist(t1[-whxx])))## selecting only those patients for which both events have valid values
#         t2 <- as.numeric(unname(unlist(t2[-whxx])))
#         c1 <- c(c1, t1)
#         c2 <- c(c2, t2)
#         p1 <- c(p1, rep(as1[j], length(t1)))
#         p2 <- c(p2, rep(as2[j], length(t2)))
#         cor1 <- c(cor1, rep(temp$COR1[j], length(t1)))
#         cor2 <- c(cor2, rep(temp$COR2[j], length(t2)))
#         tcancer <- c(tcancer, rep(all_cancer[k], length(t1)))
#     }
# }

# pdata <- data.frame(CANCER=tcancer, AS1=p1, AS2=p2, COR1=cor1, COR2=cor2, psi1=c1, psi2=c2)

# for(k in 1:length(all_cancer)){

#     # tempdata <- pdata[pdata$tcancer == all_cancer[k], ]
#     wh <- which(tcancer == all_cancer[k])

#     tempdata <- data.frame(CANCER=tcancer[wh], AS1=p1[wh], AS2=p2[wh], COR1=cor1[wh], COR2=cor2[wh], psi1=c1[wh], psi2=c2[wh])
#     tempdata$ID <- paste0(tempdata$AS1,'_',tempdata$AS2)

#     flag <- rep('1st Quadrant', length(tempdata[[1]]))
#     wh1 <- which(tempdata$COR1 > 0)
#     wh2 <- which(tempdata$COR2 < 0)
#     wha <- intersect(wh1, wh2)
#     flag[wha] <- '4th Quadrant'
#     wh1 <- which(tempdata$COR1 < 0)
#     wh2 <- which(tempdata$COR2 > 0)
#     wha <- intersect(wh1, wh2)
#     flag[wha] <- '2nd Quadrant'
#     wh1 <- which(tempdata$COR1 < 0)
#     wh2 <- which(tempdata$COR2 < 0)
#     wha <- intersect(wh1, wh2)
#     flag[wha] <- '3rd Quadrant'
#     tempdata$FLAG <- flag

#     tempdata1 <- tempdata[,c(8,9,2,4,6)]
#     colnames(tempdata1) <- c('ID','FLAG','EVENT','COR','PSI')
#     tempdata1$EVENT <- rep('X', length(tempdata1[[1]]))
#     tempdata2 <- tempdata[,c(8,9,3,5,7)]
#     colnames(tempdata2) <- c('ID','FLAG','EVENT','COR','PSI')
#     tempdata2$EVENT <- rep('Y', length(tempdata2[[1]]))
#     tempdataz <- rbind(tempdata1, tempdata2)

#     q1 <- unique(tempdataz$ID[which(tempdataz$FLAG == '1st Quadrant')])
#     q2 <- unique(tempdataz$ID[which(tempdataz$FLAG == '2nd Quadrant')])
#     q3 <- unique(tempdataz$ID[which(tempdataz$FLAG == '3rd Quadrant')])
#     q4 <- unique(tempdataz$ID[which(tempdataz$FLAG == '4th Quadrant')])
#     xlabss <- c(q1,q2,q3,q4)

#     ## for each event pair compute the median of patient-specific PSI vallues
#     md1 <- c()
#     md2 <- c()
#     for(j in 1:length(xlabss)){
#         temq <- tempdataz[tempdataz$ID == xlabss[j], ]
#         md1 <- c(md1, median(temq[temq$EVENT == 'X',]$PSI))
#         md2 <- c(md2, median(temq[temq$EVENT == 'Y',]$PSI))
#     }

#     tempdatazz <- data.frame(ID=rep(xlabss,2), PSI=c(md1,md2), EVENT=c(rep('X',length(xlabss)), rep('Y',length(xlabss))))

#     # p <- ggplot(tempdataz, aes(ID, PSI, fill=EVENT)) + 
#     # geom_boxplot()+
#     # theme(legend.text=element_text(size=12))
#     # basesize <- 12
#     # p <- p + theme_bw(base_size = basesize * 0.8) +
#     # scale_x_discrete(name="Event pair", labels=xlabss) +
#     # scale_y_continuous(name="Percent Spliced In (PSI)") +
#     # # scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
#     # geom_vline(xintercept=length(q1)+0.5, color='red')+
#     # theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
#     # axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
#     # strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
#     # guides(color=guide_legend(title="Cancer type",ncol=3))
#     # ggsave(p,filename=paste0(save_dir,"/",all_cancer[k],'_pairedPSI.png'),width=12, height=5, dpi=400)

#     p <- ggplot(tempdatazz, aes(ID, PSI, color=EVENT, group=EVENT)) + 
#     geom_point()+geom_line()+
#     theme(legend.text=element_text(size=12))
#     basesize <- 12
#     p <- p + theme_bw(base_size = basesize * 0.8) +
#     scale_x_discrete(name="Event pair", labels=xlabss) +
#     scale_y_continuous(name="Percent Spliced In (PSI)") +
#     # scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
#     geom_vline(xintercept=length(q1)+0.5, linetype=1)+
#     geom_vline(xintercept=length(q1)+length(q2)+0.5, linetype=1)+
#     geom_vline(xintercept=length(q1)+length(q2)+length(q3)+0.5, linetype=1)+
#     theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
#     axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
#     strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
#     guides(color=guide_legend(title="Cancer type",ncol=3))
#     ggsave(p,filename=paste0(save_dir,"/",all_cancer[k],'_pairedPSI.png'),width=12, height=5, dpi=400)


# }






##--- work from here ----------------

##---- take an example of events of the same gene with opposite correlation with gene expression patterns ------------------------------------
tempxx <- ttdata[ttdata$CANCER == 'BRCA',]
wh1 <- which(tempxx$COR1 < 0)
wh2 <- which(tempxx$COR2 > 0)
wha <- intersect(wh1, wh2)
wh1 <- which(tempxx$COR1 > 0)
wh2 <- which(tempxx$COR2 < 0)
whb <- intersect(wh1, wh2)
wh <- union(wha, whb)
tempyy <- tempxx[wh, ]

tempxx[order(tempxx$AS1),]





##-----------------------------------------------------------------------------------------------

##--- top 5 from each cancer ---------
whv <- c()
for(k in 1:length(all_cancer)){

    temp <- pdatay[pdatay$cancer == all_cancer[k], ]
    temp <- temp[order(-abs(temp$Spearman_value)), ]
    temp <- temp[1:5,]
    wh1 <- which(pdatay$Event %in% temp$Event)
    wh2 <- which(pdatay$cancer == all_cancer[k])

    whv <- c(whv, intersect(wh1, wh2))
}

pdataz <- pdatay[whv, ]



