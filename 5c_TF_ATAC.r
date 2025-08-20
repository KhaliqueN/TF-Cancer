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

# paired_sam <- data.table::fread('../data/cancer_paired_samples.txt')
# paired_sam <- paired_sam[-4]

input_dir <- '../data/PSI_data'
fdr <- 0.05
sample_size <- 10
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
    wha <- which(temp$FDR < fdr)
    whb <- which(abs(temp$MEAN_NORMAL-temp$MEAN_CANCER) > fdr)
    wh <- intersect(wha, whb)
    tempx <- temp[wh, ]
    
    tempy <- tempx[tempx$symbol %in% tfs$Gene_Symbol, ]
    temp_as <- tempy$as_id

    ##-- look which cancer samples (not necessarily paired) are also present for ATAC-seq
    temp_psi <- data.table::fread(all_files_raw[k], sep='\t', fill=TRUE)
    temp_psi <- as.data.frame(temp_psi)
    wh <- which(colnames(temp_psi) %in% all_samples)
    print(paste0('# of samples: ',length(wh)))

    if(length(wh) != 0){

        tempyy <- temp_psi[,c(seq(1,3),wh)]
        tempz <- tempyy[tempyy$as_id %in% temp_as, ]

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
        # tempy1 <- tempy[tempy$symbol %in% tempz$symbol,]
        # tempy1 <- tempy1[order(tempy1$as_id), ]

        tempz3 <- tempz[,c(1,2,3)]
        tempz1 <- tempz[,-c(1,2,3)]
        ## sort columns by tempq1
        tempz2 <- tempz1[match(colnames(tempq1), colnames(tempz1))] ### PSI values -----
    # print(paste0('# of events: ',length(tempz3[[1]])))
    # print(paste0('# of tfs: ',length(unique(tempz3[[1]]))))

        ##-- correlation between Footprinting depth and splicing event PSI values ----
        for(j in 1:length(tempz2[[1]])){

            temps <- as.numeric(tempz2[j,]) ## if at least 10 values are present
            tfp <- which(tempq_tf == tempz3$symbol[j])

            whx <- which(!is.na(temps))
            swhx <- sd(temps[whx])

            if(length(whx) >= sample_size & swhx != 0){
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
    cat('Cancer',k,'of',length(all_cancer),'done\n')
}

pdata <- data.frame(CANCER=tcancer, TF=ttf, ASID=tevent, SPLICE_TYPE=tstype, MOTIF=tmotif, CORR=tcor, PVAL=tpval, RND_PVAL=temp_rnd/rnd_expr)

ssum <- 0
for(k in 1:length(all_cancer)){
    temp <- pdata[pdata$CANCER == all_cancer[k], ]
    # print(length(temp[[1]]))
    temp <- temp[complete.cases(temp), ]
    # print(length(temp[[1]]))
    ssum <- ssum+length(unique(temp$ASID))
    print(paste0('# of events: ',length(unique(temp$ASID))))
    print(paste0('# of TFs: ',length(unique(temp$TF))))
}


pdata <- pdata[complete.cases(pdata), ]
TF_depth_all <- pdata
TF_depth_all$FLAG <- rep('Depth', length(TF_depth_all[[1]]))
TF_depth_all$ID <- paste0(TF_depth_all$ASID,'_',TF_depth_all$CANCER, '_', TF_depth_all$TF)

# ##--- for each ASID-cancer pair, select one motif with larger correlation --
# tempd <- data.frame(matrix(ncol=length(pdata), nrow=0))
# for(k in 1:length(all_cancer)){
#     temp <- pdata[pdata$CANCER == all_cancer[k], ]
#     temptfs <- unique(temp$ASID)
#     for(j in 1:length(temptfs)){
#         tempx <- temp[temp$ASID == temptfs[j], ]
#         tempd <- rbind(tempd, tempx[abs(tempx$CORR) == max(abs(tempx$CORR)), ])
#     }
# }

# tempd$QVAL <- p.adjust(tempd$RND_PVAL, 'fdr')
# qdata <- tempd[tempd$QVAL < fdr, ]

# ##--- plot of PSI differences and correlations with ATAC depth ----
# mdval <- rep(0, length(tempd[[1]]))

# for(k in 1:length(all_cancer)){
#     temp <- data.table::fread(all_files[k], sep='\t')
#     tempf <- tempd[tempd$CANCER == all_cancer[k], ]
#     if(nrow(tempf) != 0){
#         wh <- which(tempd$CANCER == all_cancer[k])
#         for(j in 1:length(wh)){
#             mdval[wh[j]] <- temp[temp$as_id == tempd$ASID[wh[j]], ]$MEAN_DIFF
#         }
#     }
# }

# tempd$MEAN_DIFF <- mdval
# tempd3 <- tempd



# tempd3$color <- ifelse(tempd3$QVAL < 0.05, 'a', 'b')

# basesize <- 10
# ppx <- ggplot(data = tempd3, aes(x=CORR, y=MEDIAN_DIFF, color=color)) + 
# geom_point(alpha=0.5, size=0.25)+
# scale_x_continuous()+
# scale_y_continuous()+
# # geom_hline(yintercept=0.2, color='blue', linetype='dashed')+
# # geom_hline(yintercept=-0.2, color='blue', linetype='dashed')+
# scale_color_manual(values=c('#e34a33', '#000000'))+
# # geom_vline(xintercept=0, color='blue', linetype='dashed')+
# xlab("Correlation with footprint depth")+ylab("Median of paired \nsample differences")+
# theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
#     axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
#     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"))+
# guides(color='none')
# ggsave(ppx,filename=paste0(save_dir, "/Scatter_footprint_depth.png"),width=4, height=3, dpi=400)

# ##----
# TF_depth <- tempd3[order(-abs(tempd3$CORR)), ]
# TF_depth$ID <- paste0(TF_depth$ASID,'_',TF_depth$CANCER, '_', TF_depth$TF, '_', TF_depth$MOTIF)
# # fwrite(TF_depth, paste0(save_dir, "/Scatter_footprint_depth.csv"))





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
    wha <- which(temp$FDR < fdr)
    whb <- which(abs(temp$MEAN_NORMAL-temp$MEAN_CANCER) > fdr)
    wh <- intersect(wha, whb)
    tempx <- temp[wh, ]
    tempy <- tempx[tempx$symbol %in% tfs$Gene_Symbol, ]
    temp_as <- tempy$as_id

    ##-- look which cancer samples (not necessarily paired) are also present for ATAC-seq
    temp_psi <- data.table::fread(all_files_raw[k], sep='\t', fill=TRUE)
    temp_psi <- as.data.frame(temp_psi)
    wh <- which(colnames(temp_psi) %in% all_samples)
    print(paste0('# of samples: ',length(wh)))

    if(length(wh) != 0){

        tempyy <- temp_psi[,c(seq(1,3),wh)]
        tempz <- tempyy[tempyy$as_id %in% temp_as, ]

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
        # tempy1 <- tempy[tempy$symbol %in% tempz$symbol,]
        # tempy1 <- tempy1[order(tempy1$as_id), ]

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

            if(length(whx) >= sample_size & swhx != 0){
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

    cat('Cancer',k,'of',length(all_cancer),'done\n')

}

pdata <- data.frame(CANCER=tcancer, TF=ttf, ASID=tevent, SPLICE_TYPE=tstype, MOTIF=tmotif, CORR=tcor, PVAL=tpval, RND_PVAL=temp_rnd/rnd_expr)
pdata <- pdata[complete.cases(pdata), ]
TF_flank_all <- pdata
TF_flank_all$FLAG <- rep('Flank', length(TF_flank_all[[1]]))
TF_flank_all$ID <- paste0(TF_flank_all$ASID,'_',TF_flank_all$CANCER, '_', TF_flank_all$TF)

##--- combine the depth and flanking -----
TF_depthx <- TF_depth_all#[order(TF_depth_all$ID), ]
TF_flankx <- TF_flank_all#[order(TF_flank_all$ID), ]

TF_FP_all <- rbind(TF_depthx, TF_flankx)
data.table::fwrite(TF_FP_all, paste0(save_dir,'/Footprinting_correlation.txt'),sep='\t',row.names=FALSE, quote=FALSE)

TF_FP_allx <- TF_FP_all[,-c(7,10)]
colnames(TF_FP_allx) <- c('CANCER','TF','AS_ID','AS_TYPE','DNA_MOTIF','SPEARMAN_COR','P_VALUE','FOOTPRINT','FDR')
data.table::fwrite(TF_FP_allx, paste0(save_dir,'/Footprinting_correlation.csv'),sep='\t',row.names=FALSE, quote=FALSE)

TF_FP_all <- data.table::fread(paste0(save_dir,'/Footprinting_correlation.txt'),sep='\t')
##---------------------------------------

TF_FP_all$FDR <- p.adjust(TF_FP_all$RND_PVAL, 'fdr')
# TF_FP_all_sig <- TF_FP_all[TF_FP_all$FDR < fdr, ]


##--- for each ASID-cancer pair, select one motif with larger correlation ----
TF_FP_all_dmotif <- data.frame(matrix(ncol=length(TF_FP_all), nrow=0))
TF_FP_all_fmotif <- data.frame(matrix(ncol=length(TF_FP_all), nrow=0))

for(k in 1:length(all_cancer)){

    tempc <- TF_FP_all[TF_FP_all$CANCER == all_cancer[k], ]
    if(nrow(tempc) != 0){
        tempd <- tempc[tempc$FLAG == 'Depth', ]
        temptfs <- unique(tempd$ASID)
        for(j in 1:length(temptfs)){
            tempx <- tempd[tempd$ASID == temptfs[j], ]
            TF_FP_all_dmotif <- rbind(TF_FP_all_dmotif, tempx[abs(tempx$CORR) == max(abs(tempx$CORR)), ][1,])
        }

        tempd <- tempc[tempc$FLAG == 'Flank', ]
        temptfs <- unique(tempd$ASID)
        for(j in 1:length(temptfs)){
            tempx <- tempd[tempd$ASID == temptfs[j], ]
            TF_FP_all_fmotif <- rbind(TF_FP_all_fmotif, tempx[abs(tempx$CORR) == max(abs(tempx$CORR)), ][1,])
        }
    }
}


##---- ASID-cancer pairs common in depth and flank ----------------
cmas <- intersect(TF_FP_all_dmotif$ID, TF_FP_all_fmotif$ID)


TF_FP_all_dmotif <- TF_FP_all_dmotif[TF_FP_all_dmotif$ID %in% cmas, ]
TF_FP_all_fmotif <- TF_FP_all_fmotif[TF_FP_all_fmotif$ID %in% cmas, ]

TF_depthy <- TF_FP_all_dmotif[order(TF_FP_all_dmotif$ID), ]
TF_flanky <- TF_FP_all_fmotif[order(TF_FP_all_fmotif$ID), ]


##---- Add delta PSI values -----------
mdval <- rep(0, length(TF_depthy[[1]]))
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    tempf <- TF_depthy[TF_depthy$CANCER == all_cancer[k], ]
    if(nrow(tempf) != 0){
        wh <- which(TF_depthy$CANCER == all_cancer[k])
        for(j in 1:length(wh)){
            mdval[wh[j]] <- temp[temp$as_id == TF_depthy$ASID[wh[j]], ]$MEAN_DIFF
        }
    }
}

TF_depthy$MEAN_DIFF <- mdval
TF_flanky$MEAN_DIFF <- mdval


##---- plot ---------
tempd3 <- TF_depthy
tempd3$flank_FDR <- TF_flanky$FDR
tempd3$flank_CORR <- TF_flanky$CORR

wh1 <- which(tempd3$FDR < 0.05)
wh2 <- which(tempd3$flank_FDR < 0.05)
wh3 <- intersect(wh1, wh2)
wh12 <- setdiff(wh1, wh2)
wh21 <- setdiff(wh2, wh1)

cormax <- c()
tsig <- rep('Neither', length(tempd3[[1]]))
tsig[wh12] <- 'Footprinting depth'
tsig[wh21] <- 'Flanking accessibility'
tsig[wh3] <- 'Both'

for(k in 1:length(tempd3[[1]])){
    tax <- abs(tempd3$flank_CORR[k])
    gax <- abs(tempd3$CORR[k])
    if(tax > gax){
        cormax <- c(cormax, tempd3$flank_CORR[k])
    }else{
        cormax <- c(cormax, tempd3$CORR[k])
    }
}
tempd3$CORR_MAX <- cormax
tempd3$SIG <- tsig
tempd3$IDD <- paste0(tempd3$TF,'_',tempd3$ASID,'_',tempd3$SPLICE_TYPE, '\n (', tempd3$CANCER,', \u0394PSI: ',signif(tempd3$MEAN_DIFF,2),')')
tomark <- tempd3[tempd3$SIG != 'Neither', ]$IDD ## 152 PTSE cancertype pairs are significant
# unique(tempd3[tempd3$SIG != 'Neither', ]$ASID) ## 152 PTSE cancertype pairs are significant

# tolabel1 <- subset(tempd3, IDD %in% tomark & abs(MEAN_DIFF) > 0.2)$ID
# tolabel2 <- subset(tempd3, IDD %in% tomark & abs(CORR) > 0.5)$ID
# tolabel3 <- subset(tempd3, IDD %in% tomark & abs(flank_CORR) > 0.5)$ID

# tolabela <- intersect(tolabel1, tolabel2)
# tolabelb <- intersect(tolabel1, tolabel3)
# tolabel <- union(tolabela, tolabelb)
tolabel <- subset(tempd3, IDD %in% tomark & abs(MEAN_DIFF) > 0.4)$ID

whl <- which(tempd3$ID %in% tolabel)
tolabelx <- tempd3$IDD[which(tempd3$ID %in% tolabel)]
tempd3$PL <- ""
tempd3$PL[whl] <- tolabelx
tempd3 <- tempd3 %>% arrange(desc(SIG))

basesize <- 10
ppx <- ggplot(data = tempd3, aes(x=flank_CORR, y=CORR, color=SIG, label=PL)) + 
geom_point(alpha=0.8, size=1)+
scale_x_continuous(limits=c(-1.4,1.4))+
scale_y_continuous(limits=c(-1.4,1.4))+
scale_color_manual(values=c('#1b9e77', '#d95f02', '#7570b3', '#d9d9d9'))+
xlab("Correlation with flanking accessibility")+ylab("Correlation with footprinting depth")+
geom_text_repel(family = "Poppins",
    max.overlaps=Inf,
                      size = 3,
                      color='black',
                      arrow = arrow(length = unit(0.010, "npc")),
                      min.segment.length = 0) +
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(color=guide_legend(title="Significance",ncol=1,override.aes = list(size = 2)))
ggsave(ppx,filename=paste0(save_dir, "/Footprint_depth_flank_sig.png"),width=5.5, height=3.5, dpi=400)


###--------------------------------------------------------------------------
