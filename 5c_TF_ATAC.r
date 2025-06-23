##############################################################################################
# Purpose: Overlap of significant TF splicing events with cancer type-related footprinting
## Do it for all perturbed events and not just the ones with DBD perturbations
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(GenomicDataCommons)

save_dir <- '../results_new/Footprinting'
if(dir.exists(save_dir)){
    unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

##-- TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')

# dbd_purt <- data.table::fread('../data/Events_perturbing_DBD.txt')
# ed_purt <- data.table::fread('../data/Events_perturbing_ED.txt')
##----------------------------------------------------------------------

paired_sam <- data.table::fread('../data/cancer_paired_samples.txt')
paired_sam <- paired_sam[-4]

input_dir <- '../data/PSI_data'
fdr <- 0.05
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))
all_files <- all_files[-4]
all_files_raw <- all_files_raw[-4]
all_cancer <- substr(basename(all_files), 1,4)



## Footprinting depth data ---------------
## From this study: The chromatin accessibility landscape of primary human cancers ----
tmas <- as.data.frame(readxl::read_excel('../data/aav1898_data_s6.xlsx', 1))
tmasx <- tmas[17:length(tmas[[1]]),]
colnames(tmasx) <- tmasx[1,]
tmasx <- tmasx[-1,]

# tmas2 <- as.data.frame(readxl::read_excel('../data/aav1898_data_s6.xlsx', 3))
# tmasx2 <- tmas2[17:length(tmas2[[1]]),]
# ## filter the first file to only retain these motifs
# tmasx <- tmasx[tmasx$`CIS-BP_ID` %in% tmasx2[[2]], ]

all_samples <- colnames(tmasx)[-c(1,2,3,4)]
all_samples <- gsub('-','_', all_samples)
rnd_expr <- 100

tcancer <- c()
tcor <- c()
tmotif <- c()
tpval <- c()
ttf <- c()
tevent <- c()
tstype <- c()
temp_rnd <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    # temp_as <- dbd_purt[dbd_purt$CANCER == all_cancer[k], ]$AS
    # tempy <- tempx[tempx$as_id %in% temp_as, ]
    tempy <- tempx[tempx$symbol %in% tfs$Gene_Symbol, ]
    temp_as <- tempy$as_id

    ##-- look which cancer samples (not necessarily paired) are also present for ATAC-seq
    temp_psi <- data.table::fread(all_files_raw[k], sep='\t', fill=TRUE)
    temp_psi <- as.data.frame(temp_psi)
    wh <- which(colnames(temp_psi) %in% all_samples)

    if(length(wh) != 0){

        tempyy <- temp_psi[,c(seq(1,3),wh)]
        tempz <- tempyy[tempyy$as_id %in% temp_as, ]
        # tempz <- tempz[order(tempz$as_id), ]

        # tempz$MEDIAN_DIFF <- tempy$MEDIAN_DIFF
        # tempz$FDR <- tempy$FDR

        ##-- ATAC-seq depth ---
        wh <- which(all_samples %in% intersect(colnames(temp_psi), all_samples))+4
        tempq <- tmasx[tmasx$TF_Name %in% tempz$symbol, c(2,3,wh)]
        tempq_tf <- tempq[[2]]
        tempq_motif <- tempq[[1]]
        tempq1 <- tempq[,-c(1,2)]
        colnames(tempq1) <- gsub('-','_',colnames(tempq1)) ### Footprinting depth -----

        ##-- filter tempz based on available atac-seq data ---
        tempz <- tempz[tempz$symbol %in% tempq$TF_Name, ]
        tempz <- tempz[order(tempz$as_id), ]
        tempy1 <- tempy[tempy$symbol %in% tempz$symbol,]
        tempy1 <- tempy1[order(tempy1$as_id), ]

        tempz3 <- tempz[,c(1,2,3)]
        tempz1 <- tempz[,-c(1,2,3)]
        ## sort columns by tempq1
        tempz2 <- tempz1[match(colnames(tempq1), colnames(tempz1))] ### PSI values -----

        ##-- correlation between Footprinting depth and splicing event PSI values ----
        for(j in 1:length(tempz2[[1]])){

            temps <- as.numeric(tempz2[j,]) ## if at least 5 values are present
            tfp <- which(tempq_tf == tempz3$symbol[j])

            whx <- which(!is.na(temps))
            swhx <- sd(temps[whx])

            if(length(whx) > 5 & swhx != 0){
                for(i in 1:length(tfp)){
                    tempqs <- as.numeric(tempq1[tfp[i], ])
                    temp_cor <- cor.test(tempqs, temps, method='spearman', use="complete.obs")
                    t_rnd_cors <- c()
                    for(ii in 1:rnd_expr){ ##--- randomized experiment
                        temp_cor_rnd <- cor.test(x=tempqs, y=sample(temps), method = 'spearman', use="complete.obs")
                        t_rnd_cors <- c(t_rnd_cors, as.numeric(temp_cor_rnd$estimate))
                    }

                    if(temp_cor$estimate > 0){
                        temp_rnd <- c(temp_rnd, length(which(t_rnd_cors > as.numeric(temp_cor$estimate))))
                    }else{
                        temp_rnd <- c(temp_rnd, length(which(t_rnd_cors < as.numeric(temp_cor$estimate))))
                    }

                    tcancer <- c(tcancer, all_cancer[k])
                    tcor <- c(tcor, temp_cor$estimate)
                    tmotif <- c(tmotif, tempq_motif[tfp[i]])
                    tpval <- c(tpval, temp_cor$p.value)
                    ttf <- c(ttf, tempq_tf[tfp[i]])
                    tevent <- c(tevent, tempz3$as_id[j])
                    tstype <- c(tstype, tempz3$splice_type[j])
                }
            }
        }
    }
}

# pdata <- data.frame(CANCER=tcancer, PSI=tpsi, FPD=tfpd)
pdata <- data.frame(CANCER=tcancer, TF=ttf, ASID=tevent, SPLICE_TYPE=tstype, MOTIF=tmotif, CORR=tcor, PVAL=tpval, RND_PVAL=temp_rnd/rnd_expr)
pdata <- pdata[complete.cases(pdata), ]
# pdata$ID <- paste0(pdata$ASID,'_',pdata$CANCER, '_', pdata$TF, '_', pdata$MOTIF)

# pdata$QVAL <- p.adjust(pdata$RND_PVAL, 'fdr')
# qdata <- pdata[pdata$QVAL < 0.05, ]

##--- for each ASID-cancer pair, select one motif with larger correlation --
tempd <- data.frame(matrix(ncol=length(pdata), nrow=0))
for(k in 1:length(all_cancer)){
    temp <- pdata[pdata$CANCER == all_cancer[k], ]
    temptfs <- unique(temp$ASID)
    for(j in 1:length(temptfs)){
        tempx <- temp[temp$ASID == temptfs[j], ]
        tempd <- rbind(tempd, tempx[abs(tempx$CORR) == max(abs(tempx$CORR)), ])
    }
}

tempd$QVAL <- p.adjust(tempd$RND_PVAL, 'fdr')
qdata <- tempd[tempd$QVAL < 0.05, ]

##--- plot of PSI differences and correlations with ATAC depth ----
mdval <- rep(0, length(tempd[[1]]))

for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    tempf <- tempd[tempd$CANCER == all_cancer[k], ]
    if(nrow(tempf) != 0){
        wh <- which(tempd$CANCER == all_cancer[k])
        for(j in 1:length(wh)){
            mdval[wh[j]] <- temp[temp$as_id == tempd$ASID[wh[j]], ]$MEDIAN_DIFF
        }
    }
}

tempd$MEDIAN_DIFF <- mdval

tempd3 <- tempd
tempd3$color <- ifelse(tempd3$QVAL < 0.05, 'a', 'b')

basesize <- 10
ppx <- ggplot(data = tempd3, aes(x=CORR, y=MEDIAN_DIFF, color=color)) + 
geom_point(alpha=0.5, size=0.25)+
scale_x_continuous()+
scale_y_continuous()+
# geom_hline(yintercept=0.2, color='blue', linetype='dashed')+
# geom_hline(yintercept=-0.2, color='blue', linetype='dashed')+
scale_color_manual(values=c('#e34a33', '#000000'))+
# geom_vline(xintercept=0, color='blue', linetype='dashed')+
xlab("Correlation with footprint depth")+ylab("Median of paired \nsample differences")+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(color='none')
ggsave(ppx,filename=paste0(save_dir, "/Scatter_footprint_depth.png"),width=4, height=3, dpi=400)

##----
TF_depth <- tempd3[order(-abs(tempd3$MEDIAN_DIFF)), ]
TF_depth$ID <- paste0(TF_depth$ASID,'_',TF_depth$CANCER, '_', TF_depth$TF, '_', TF_depth$MOTIF)
fwrite(TF_depth, paste0(save_dir, "/Scatter_footprint_depth.csv"))





## Footprinting flanking data ---------------
## From this study: The chromatin accessibility landscape of primary human cancers ----
tmas <- as.data.frame(readxl::read_excel('../data/aav1898_data_s6.xlsx', 2))
tmasx <- tmas[17:length(tmas[[1]]),]
colnames(tmasx) <- tmasx[1,]
tmasx <- tmasx[-1,]

all_samples <- colnames(tmasx)[-c(1,2,3,4)]
all_samples <- gsub('-','_', all_samples)
rnd_expr <- 100

tcancer <- c()
tcor <- c()
tmotif <- c()
tpval <- c()
ttf <- c()
tevent <- c()
tstype <- c()
temp_rnd <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    # temp_as <- dbd_purt[dbd_purt$CANCER == all_cancer[k], ]$AS
    # tempy <- tempx[tempx$as_id %in% temp_as, ]
    tempy <- tempx[tempx$symbol %in% tfs$Gene_Symbol, ]
    temp_as <- tempy$as_id

    ##-- look which cancer samples (not necessarily paired) are also present for ATAC-seq
    temp_psi <- data.table::fread(all_files_raw[k], sep='\t', fill=TRUE)
    temp_psi <- as.data.frame(temp_psi)
    wh <- which(colnames(temp_psi) %in% all_samples)

    if(length(wh) != 0){

        tempyy <- temp_psi[,c(seq(1,3),wh)]
        tempz <- tempyy[tempyy$as_id %in% temp_as, ]
        # tempz <- tempz[order(tempz$as_id), ]

        # tempz$MEDIAN_DIFF <- tempy$MEDIAN_DIFF
        # tempz$FDR <- tempy$FDR

        ##-- ATAC-seq depth ---
        wh <- which(all_samples %in% intersect(colnames(temp_psi), all_samples))+4
        tempq <- tmasx[tmasx$TF_Name %in% tempz$symbol, c(2,3,wh)]
        tempq_tf <- tempq[[2]]
        tempq_motif <- tempq[[1]]
        tempq1 <- tempq[,-c(1,2)]
        colnames(tempq1) <- gsub('-','_',colnames(tempq1)) ### Footprinting depth -----

        ##-- filter tempz based on available atac-seq data ---
        tempz <- tempz[tempz$symbol %in% tempq$TF_Name, ]
        tempz <- tempz[order(tempz$as_id), ]
        tempy1 <- tempy[tempy$symbol %in% tempz$symbol,]
        tempy1 <- tempy1[order(tempy1$as_id), ]

        tempz3 <- tempz[,c(1,2,3)]
        tempz1 <- tempz[,-c(1,2,3)]
        ## sort columns by tempq1
        tempz2 <- tempz1[match(colnames(tempq1), colnames(tempz1))] ### PSI values -----

        ##-- correlation between Footprinting depth and splicing event PSI values ----
        for(j in 1:length(tempz2[[1]])){

            temps <- as.numeric(tempz2[j,]) ## if at least 5 values are present
            tfp <- which(tempq_tf == tempz3$symbol[j])

            whx <- which(!is.na(temps))
            swhx <- sd(temps[whx])

            if(length(whx) > 5 & swhx != 0){
                for(i in 1:length(tfp)){
                    tempqs <- as.numeric(tempq1[tfp[i], ])
                    temp_cor <- cor.test(tempqs, temps, method='spearman', use="complete.obs")
                    t_rnd_cors <- c()
                    for(ii in 1:rnd_expr){ ##--- randomized experiment
                        temp_cor_rnd <- cor.test(x=tempqs, y=sample(temps), method = 'spearman', use="complete.obs")
                        t_rnd_cors <- c(t_rnd_cors, as.numeric(temp_cor_rnd$estimate))
                    }

                    if(temp_cor$estimate > 0){
                        temp_rnd <- c(temp_rnd, length(which(t_rnd_cors > as.numeric(temp_cor$estimate))))
                    }else{
                        temp_rnd <- c(temp_rnd, length(which(t_rnd_cors < as.numeric(temp_cor$estimate))))
                    }

                    tcancer <- c(tcancer, all_cancer[k])
                    tcor <- c(tcor, temp_cor$estimate)
                    tmotif <- c(tmotif, tempq_motif[tfp[i]])
                    tpval <- c(tpval, temp_cor$p.value)
                    ttf <- c(ttf, tempq_tf[tfp[i]])
                    tevent <- c(tevent, tempz3$as_id[j])
                    tstype <- c(tstype, tempz3$splice_type[j])
                }
            }
        }
    }
}


# pdata <- data.frame(CANCER=tcancer, PSI=tpsi, FPD=tfpd)
pdata <- data.frame(CANCER=tcancer, TF=ttf, ASID=tevent, SPLICE_TYPE=tstype, MOTIF=tmotif, CORR=tcor, PVAL=tpval, RND_PVAL=temp_rnd/rnd_expr)
pdata <- pdata[complete.cases(pdata), ]
pdata$ID <- paste0(pdata$ASID,'_',pdata$CANCER, '_', pdata$TF, '_', pdata$MOTIF)

##--- for each TF-cancer pair, select one motif with larger correlation --
tempd <- pdata[pdata$ID %in% TF_depth$ID, ]

tempd$QVAL <- p.adjust(tempd$RND_PVAL, 'fdr')
qdata <- tempd[tempd$QVAL < 0.05, ]


##--
TF_depthx <- TF_depth[order(TF_depth$ID), ]
tempd <- tempd[order(tempd$ID), ]

tempd$DEPTH <- TF_depthx$CORR
tempd$DEPTH_QVAL <- TF_depthx$QVAL
tempd$MEDIAN_DIFF <- TF_depthx$MEDIAN_DIFF


tempd3 <- tempd
wh1 <- which(tempd3$DEPTH_QVAL < 0.05)
wh2 <- which(tempd3$QVAL < 0.05)
wh3 <- intersect(wh1, wh2)
wh12 <- setdiff(wh1, wh2)
wh21 <- setdiff(wh2, wh1)

cormax <- c()
tsig <- rep('Neither', length(tempd3[[1]]))
tsig[wh12] <- 'Depth'
tsig[wh21] <- 'Flank'
tsig[wh3] <- 'Both'

for(k in 1:length(tempd3[[1]])){
    tax <- abs(tempd3$CORR[k])
    gax <- abs(tempd3$DEPTH[k])
    if(tax > gax){
        cormax <- c(cormax, tempd3$CORR[k])
    }else{
        cormax <- c(cormax, tempd3$DEPTH[k])
    }
}
tempd3$CORR_MAX <- cormax
tempd3$SIG <- tsig

tempd4 <- tempd3[union(wh1, wh2), ]
# tempd3$CORR_MAX <- mapply(function(x,y) max(abs(x), abs(y)), tempd3$CORR, tempd3$DEPTH)

basesize <- 10
ppx <- ggplot(data = tempd3, aes(x=CORR_MAX, y=MEDIAN_DIFF)) + 
geom_point(alpha=0.6)+
scale_x_continuous()+
scale_y_continuous()+
geom_hline(yintercept=0.2, color='red', linetype='dashed')+
geom_hline(yintercept=-0.2, color='red', linetype='dashed')+
# scale_color_manual(values=c('#e34a33', '#000000'))+
# geom_vline(xintercept=0, color='blue', linetype='dashed')+
xlab("Correlation with footprint depth/flank")+ylab("Median of paired \nsample differences")+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(color='none')
ggsave(ppx,filename=paste0(save_dir, "/Scatter_footprint_depth_flank.png"),width=4, height=3, dpi=400)

ppx <- ggplot(data = tempd4, aes(x=CORR_MAX, y=MEDIAN_DIFF)) + 
geom_point(alpha=0.6)+
scale_x_continuous()+
scale_y_continuous()+
geom_hline(yintercept=0.2, color='red', linetype='dashed')+
geom_hline(yintercept=-0.2, color='red', linetype='dashed')+
# scale_color_manual(values=c('#e34a33', '#000000'))+
# geom_vline(xintercept=0, color='blue', linetype='dashed')+
xlab("Correlation with footprint depth/flank")+ylab("Median of paired \nsample differences")+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(color='none')
ggsave(ppx,filename=paste0(save_dir, "/Scatter_footprint_depth_flank_sig.png"),width=4, height=3, dpi=400)


tempd3 <- tempd3 %>% arrange(desc(SIG))
tempd3$IDD <- paste0(tempd3$TF,'_',tempd3$ASID,'_',tempd3$SPLICE_TYPE, ' (',signif(tempd3$MEDIAN_DIFF,2),')')
# tolabel1 <- subset(tempd3, abs(CORR) > 0.75 & QVAL < 0.05 & abs(MEDIAN_DIFF) > 0.1)$ID
# tolabel2 <- subset(tempd3, abs(DEPTH) > 0.75 & DEPTH_QVAL < 0.05 & abs(MEDIAN_DIFF) > 0.1)$ID
tolabel1 <- subset(tempd3, QVAL < 0.05 & abs(MEDIAN_DIFF) > 0.3)$ID
tolabel2 <- subset(tempd3, DEPTH_QVAL < 0.05 & abs(MEDIAN_DIFF) > 0.3)$ID
tolabel <- union(tolabel1, tolabel2)
whl <- which(tempd3$ID %in% tolabel)
tolabelx <- tempd3$IDD[which(tempd3$ID %in% tolabel)]
tempd3$PL <- ""
tempd3$PL[whl] <- tolabelx

ppx <- ggplot(data = tempd3, aes(x=CORR, y=DEPTH, color=SIG, label=PL)) + 
geom_point(alpha=0.8, size=0.5)+
scale_x_continuous(limits=c(-1.5,1.5))+
scale_y_continuous(limits=c(-1.5,1.5))+
scale_color_manual(values=c('#1b9e77', '#d95f02', '#7570b3', '#d9d9d9'))+
xlab("Correlation with footprint flank")+ylab("Correlation with footprint depth")+
geom_text_repel( family = "Poppins",
    size = 2,
    min.segment.length = 0, 
    seed = 42, 
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .15,
    nudge_y = .5,
    color = "black") +
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(color=guide_legend(title="Significance",ncol=1,override.aes = list(size = 2)))
ggsave(ppx,filename=paste0(save_dir, "/Footprint_depth_flank_sig.png"),width=5, height=3.5, dpi=400)


##-----------
tempd5 <- tempd3[order(-abs(tempd3$MEDIAN_DIFF)), ]
TF_depth_flank <- tempd5
fwrite(TF_depth_flank, paste0(save_dir, "/Scatter_footprint_depth_flank.csv"))






##------ TF depth/flanks not correlated with gene expression but with splicing behaviour -----
tmas <- as.data.frame(readxl::read_excel('../data/aav1898_data_s6.xlsx', 3))
tmasx <- tmas[17:length(tmas[[1]]),]
colnames(tmasx) <- tmasx[1,]
tmasx <- tmasx[-1,]
cor_tf <- tmasx$TF_Name[union(which(tmasx$depth_FDR < 0.05), which(tmasx$flank_FDR < 0.05))]
# ncor_tf <- tmasx$TF_Name[union(which(tmasx$depth_FDR > 0.5), which(tmasx$flank_FDR > 0.5))]
ncor_tf <- tmasx$TF_Name[intersect(which(abs(as.numeric(tmasx$depth_correlation)) < 0.1), which(abs(as.numeric(tmasx$flank_correlation)) < 0.1))]

tempd5_cor <- tempd5[tempd5$TF %in% cor_tf, ]
wh1 <- which(tempd5_cor$QVAL < 0.05)
wh2 <- which(tempd5_cor$DEPTH_QVAL < 0.05)
wh <- union(wh1, wh2)
tempd5_cor <- tempd5_cor[wh,]
tempd5_cor <- tempd5_cor[order(-abs(tempd5_cor$MEDIAN_DIFF)),]
tempd5_cor$FLAG <- 'RNA-seq and Splicing'

tempd5_ncor <- tempd5[tempd5$TF %in% ncor_tf, ]
wh1 <- which(tempd5_ncor$QVAL < 0.05)
wh2 <- which(tempd5_ncor$DEPTH_QVAL < 0.05)
wh <- union(wh1, wh2)
tempd5_ncor <- tempd5_ncor[wh,]
tempd5_ncor <- tempd5_ncor[order(-abs(tempd5_ncor$MEDIAN_DIFF)),]
tempd5_ncor$FLAG <- 'Splicing'

tempd6 <- rbind(tempd5_cor, tempd5_ncor)

tempd6$IDD <- paste0(tempd6$TF,'_',tempd6$ASID,'_',tempd6$SPLICE_TYPE, ' (',signif(tempd6$MEDIAN_DIFF,2),')')
tolabel1 <- subset(tempd6, QVAL < 0.05 & abs(MEDIAN_DIFF) > 0.3)$ID
tolabel2 <- subset(tempd6, DEPTH_QVAL < 0.05 & abs(MEDIAN_DIFF) > 0.3)$ID
tolabel <- union(tolabel1, tolabel2)
whl <- which(tempd3$ID %in% tolabel)
tolabelx <- tempd6$IDD[which(tempd6$ID %in% tolabel)]
tempd6$PL <- ""
tempd6$PL[whl] <- tolabelx

ppx <- ggplot(data = tempd6, aes(x=CORR, y=DEPTH, color=FLAG, label=PL)) + 
geom_point(alpha=0.8, size=0.5)+
scale_x_continuous(limits=c(-1.5,1.5))+
scale_y_continuous(limits=c(-1.5,1.5))+
scale_color_manual(values=c('#1b9e77', '#d95f02', '#7570b3', '#d9d9d9'))+
xlab("Correlation with footprint flank")+ylab("Correlation with footprint depth")+
geom_text_repel( family = "Poppins",
    size = 2,
    min.segment.length = 0, 
    seed = 42, 
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .15,
    nudge_y = .5,
    color = "black") +
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(color=guide_legend(title="Significance",ncol=1,override.aes = list(size = 2)))
ggsave(ppx,filename=paste0(save_dir, "/Footprint_depth_flank_diff.png"),width=5, height=3.5, dpi=400)



# for(k in 1:length(all_cancer)){
#     pdatax <- pdata[pdata$CANCER == all_cancer[k], ]
#     if(nrow(pdatax) != 0){
#         pdatax <- pdatax[complete.cases(pdatax),]
#         ##--- correlation plot ----
#         p <- ggplot(pdatax, aes(PSI, FPD)) + 
#         geom_point()+
#         theme(legend.text=element_text(size=12))
#         basesize <- 12
#         p <- p + theme_bw(base_size = basesize * 0.8) +
#         scale_x_continuous(name="PSI") + 
#         scale_y_continuous(name="Footprint depth") +
#         # geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
#         # geom_text(aes(label=value), position=position_stack(vjust=0.5), size=3)+
#         theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
#         axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
#         panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#         strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
#         guides(fill='none')
#         ggsave(p,filename=paste0(save_dir,"/PSI_FPD_corr_",all_cancer[k],".png"),width=5, height=5, dpi=400)
#     }
    
# }









