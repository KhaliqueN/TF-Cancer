##############################################################################################
rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)
library(ggrepel)
library(dplyr)

save_dir <- 'results_rep/TF_splicing'

if(dir.exists(save_dir)){
	unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

# ###---------------------------------------

## TFs -------------------
tfs <- data.table::fread('data/filtered_TFs_curated.txt', sep='\t')
paired_sam <- data.table::fread('data/cancer_paired_samples.txt')
##----------------------------------------------------------------------

input_dir <- 'data/PSI_data'
fdr <- 0.05
# num_sams <- 0.5
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*PSI_paired.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))
all_cancer <- substr(basename(all_files), 1,4)

##------------------------------------------------
events_tf<- list()
tcancer <- c()
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()
RI <- c()
tumap <- list()
num_tfs <- list()
num_sig_events <- list()
num_sig_events_tfs <- list()
num_all_events <- list()
num_all_events_tfs <- list()

wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Supplementary_Table_S1.xlsx'))
##--- column name explanations --
tdatat <- data.frame(`Table fields explanation`=c('SYMBOL: HGNC gene symbol',
   'AS_ID: Alternative splicing ID provided by TCGA SpliceSeq',
   'AS_TYPE: Type of alternative splicing',
   'P_VALUE: wilcoxon signed-rank test p-value for comparison of PSI values of the concerned 
   splicing event in paired samples of the concerned cancer type',
   'FDR: Benjamini Hochberg correction of the P_VALUES for all splicing events tested',
   'MEDIAN_DIFF: Median of the differences of PSI values between paired normal and cancer samples',
     'MEDIAN_NORMAL: Median of the PSI values of normal samples',
     'MEDIAN_CANCER: Median of the PSI values of cancer samples',
    'MEAN_DIFF: Mean of the differences of PSI values between paired normal and cancer samples',
     'MEAN_NORMAL: Mean of the PSI values of normal samples',
     'MEAN_CANCER: Mean of the PSI values of cancer samples', 
     'CANCER_HIGH: Number of cancer samples in which the concerned splicing event has a higher PSI value than in normal samples',
     'NORMAL_HIGH: Number of normal samples in which the concerned splicing event has a higher PSI value than in cancer samples',
     'PAIRED_SAMPLES: Total number of paired samples of the concerned cancer type'
     ))
openxlsx::addWorksheet(wb1, sheetName = 'INDEX')
openxlsx::writeData(wb1, sheet = 'INDEX', tdatat)
openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S1.xlsx'), overwrite = T)

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    num_all_events[[k]] <- temp$as_id
    whx <- which(toupper(temp$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    temp1 <- temp[whx,]
    num_all_events_tfs[[k]] <- temp1$as_id

    wha <- which(temp$FDR < fdr)
    whb <- which(abs(temp$MEAN_NORMAL-temp$MEAN_CANCER) > fdr)
    wh <- intersect(wha, whb)
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]

    num_sig_events[[k]] <- tempx$as_id
    num_sig_events_tfs[[k]] <- tempy$as_id

    num_tfs[[k]] <- unique(tempy$symbol)

    ## save excel sheet ----
    tempz <- data.frame(tempy)
    tdatat <- tempz[, c(1,2,3, seq(length(tempz)-9, length(tempz)))]
    tdatat <- tdatat[order(-abs(tdatat$MEAN_DIFF)), ]
    tdatat$PAIRED_SAMPLES <- rep(paired_sam$PAIRED_SAMPLES[k], length(tdatat[[1]]))

    colnames(tdatat) <- c('SYMBOL', 'AS_ID', 'AS_TYPE', 'P_VALUE', 'FDR',
    'MEDIAN_DIFF', 'MEDIAN_NORMAL', 'MEDIAN_CANCER', 'MEAN_DIFF', 'MEAN_NORMAL', 
    'MEAN_CANCER', 'CANCER_HIGH', 'NORMAL_HIGH', 'PAIRED_SAMPLES' )
    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tdatat)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S1.xlsx'), overwrite = T)

    events_tf[[k]] <- unique(tempy$as_id)
    
    temp_count <- plyr::count(tempy$splice_type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RI <- c(RI, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))

    ##--- for umap ---
    whu <- which(toupper(temp$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempu <- temp[whu,]
    tumap[[k]] <- tempu$as_id

}

##--- plot the number of splicing events affecting TFs -----------------
pdatae <- data.frame(cancer=all_cancer, count1=lengths(events_tf), count2=lengths(num_tfs))#, count3=lengths(num_sig_events))
maxv <- max(pdatae$count1)
colnames(pdatae) <- c('Cancer','PTSE','TF')
pdata <- reshape2::melt(pdatae)
pdata[pdata == 0] <- NA
p <- ggplot(pdata, aes(Cancer, value, fill=variable)) + 
geom_bar(stat="identity",position="dodge")
basesize <- 8
p <- p + theme_classic(base_size = basesize) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of PTSEs or TFs", limits=c(0,maxv+100)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
scale_fill_manual("Entity",values=c('#d95f02','#1b9e77'))+
# scale_pattern_manual("Entity",values = c("crosshatch", "stripe"))+
geom_text(aes(label=value), position=position_dodge(width=0.9),hjust=0, vjust=0.5, angle=85, size=3)+
theme(axis.text.x = element_text(size = basesize, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize), axis.title=element_text(size=basesize), legend.position=c(0.85,0.86))+
guides(fill=guide_legend(title="Entity",ncol=1))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs.png"),width=3.5, height=2.8, dpi=500)

# ##------------------------------------------------------------------------------


##--- plot the number of splicing events of different types affecting TFs -----------------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=lengths(events_tf), Cancer=all_cancer)
pData <- pData[pData$Cancer != 'ESCA', ]
pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
# geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,101), breaks = seq(0, 100, by = 20))+
xlab("Cancer type")+ylab("% of PTSEs")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_classic()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"), legend.position="bottom",axis.title=element_text(size=basesize),
    legend.text=element_text(size=basesize), legend.title=element_text(size=basesize), legend.key.height= unit(0.25, 'cm'),
    legend.key.width= unit(0.25, 'cm'))+
guides(fill=guide_legend(title="Type of\nalternative\nsplicing event", ncol=1))
ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_TFs.png"),width=3, height=3.2, dpi=600)


## number of TFs --
# num_tfs_all <- unique(unlist(num_TFs))
# num_tfs_all <- unique(unlist(num_sig_events_tfs))

# ###------------------------- Correlation with number of mutated splicing factors with a mutation rate of atleast 1% -----
# ## Splicing factor mutation rate data from this study (Table S1): Somatic Mutational Landscape of Splicing Factor Genes and Their Functional Consequences across 33 Cancer Types ----
# tmas <- as.data.frame(readxl::read_excel('data/1-s2.0-S2211124718301529-mmc2.xlsx', 1))
# tmas1 <- tmas[413:length(tmas[[1]]), c(1,2,3,8,9)]
# colnames(tmas1) <- c('CANCER', 'GENE', 'MRATE', 'HS', 'LOF')
# tmas2 <- tmas1[tmas1$CANCER %in% all_cancer, ]
# tmas2[is.na(tmas2)] <- 0

# frac1 <- seq(0,50,5)
# frac2 <- seq(0,10,1)

# cor_lof <- c()
# cor_hs <- c()
# cor_mr <- c()

# for(j in 1:length(frac1)){
#     mcount_lof <- c()
#     mcount_hs <- c()
#     mcount_mr <- c()
#     for(k in 1:length(all_cancer)){
#         temp <- tmas2[tmas2$CANCER == all_cancer[k], ]
#         tempx <- temp[temp$LOF >= frac1[j], ]
#         mcount_lof <- c(mcount_lof, length(tempx[[1]]))
#         tempx <- temp[temp$HS >= frac1[j], ]
#         mcount_hs <- c(mcount_hs, length(tempx[[1]]))
#         tempx <- temp[temp$MRATE >= frac2[j], ]
#         mcount_mr <- c(mcount_mr, length(tempx[[1]]))
#     }

#     ##- correlation--
#     coral <- cor.test(x=lengths(num_sig_events), y=mcount_lof, method = 'spearman')
#     cor_lof <- c(cor_lof, coral$estimate[[1]])

#     coral <- cor.test(x=lengths(num_sig_events), y=mcount_hs, method = 'spearman')
#     cor_hs <- c(cor_hs, coral$estimate[[1]])

#     coral <- cor.test(x=lengths(num_sig_events), y=mcount_mr, method = 'spearman')
#     cor_mr <- c(cor_mr, coral$estimate[[1]])

# }


# ##--correlation plot --
# coral <- cor.test(x=pdata$Control, y=pdata$Survival, method = 'spearman')
# p <- ggplot(pdata, aes(Control, Survival, color=Cancer)) + 
# geom_point(size=2)+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_continuous(name="# of control edges") + 
# scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
# geom_text(aes(x=x1[ntype],y=y1[ntype], label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
# scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(color=guide_legend(title="Cancer type",ncol=3))
# ggsave(p,filename=paste0(save_dir,"/",net_type[ntype],"_correlation_control.png"),width=4.8, height=2.5, dpi=400)




###----- UMAP/TSNE ------------------------------------------------------------------------------------
casidx <- tumap[[1]]
for(k in 2:length(all_cancer)){
    casidx <- intersect(casidx, tumap[[k]])
} ## casidx contains TFs present in each cancer

alltfe <- events_tf[[1]]
for(k in 2:length(all_cancer)){
   alltfe <- union(alltfe, events_tf[[k]])
}

casid <- intersect(casidx, alltfe)


umap_data <- data.frame(matrix(ncol=0, nrow=0))
csamp <- c()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files_raw[k], sep='\t', fill=TRUE)
    tempx <- temp[temp$as_id %in% casid, ]
    tempx <- as.data.frame(tempx[order(tempx$as_id), ])
    tempy <- tempx[,which(colnames(tempx) %like% 'TCGA')]
    tempy <- sapply(tempy, as.numeric)
    tempy <- as.data.frame(t(tempy))
    wha <- which(rownames(tempy) %like% 'Norm')
    nflag <- rep('Cancer',length(tempy[[1]]))
    nflag[wha] <- 'Normal'
    cflag <- rep(all_cancer[k], length(tempy[[1]]))
    tempy$Sample <- nflag
    tempy$Cancer <- cflag
    tempyy <- tempy[tempy$Sample == 'Cancer', ]
    umap_data <- rbind(umap_data, tempyy)
    csamp <- c(csamp, length(tempyy[[1]]))
}

umap_data_na <- umap_data[is.na(umap_data)] <- 0
tounmap <- umap_data[, grep("V", colnames(umap_data))]

##--- TSNE ---
umap_reduction <- Rtsne::Rtsne(tounmap, check_duplicates = FALSE)
umap_reduction_df <- data.frame(V1=umap_reduction$Y[,1], V2=umap_reduction$Y[,2])
umap_reduction_df$Sample <- umap_data$Sample
umap_reduction_df$Cancer <- umap_data$Cancer

##--- tsne plot -----------
umap_reduction_df_c <- umap_reduction_df[umap_reduction_df$Sample == 'Cancer',]
p <- ggplot(umap_reduction_df_c, aes(V1, V2, color=Cancer)) + 
geom_point(alpha=0.75, size=0.5)+
theme(legend.text=element_text(size=8))
basesize <- 8
p <- p + theme_classic(base_size = basesize) +
scale_x_continuous(name="tSNE dimension 1") + 
scale_y_continuous(name="tSNE dimension 2") +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99',
    '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# geom_text(aes(label=value), position=position_stack(vjust=0.5), size=3)+
theme(axis.text.x = element_text(size = basesize*1, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize*1, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize*1), axis.title=element_text(size=basesize*1), legend.position="right", 
legend.text=element_text(size=basesize*1), legend.title=element_text(size=basesize*1))+
guides(color=guide_legend(title="Cancer type", ncol=2, override.aes = list(size = 2)))
ggsave(p,filename=paste0(save_dir,"/tSNE.png"),width=3.5, height=3.2, dpi=500)


###------ Distribution of perturbation values ----------------------------------------------------
##------------------------------------------------
tpos <- c()
tneg <- c()
tcancer <- c()
tdiff <- c()
asid <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wha <- which(temp$FDR < fdr)
    whb <- which(abs(temp$MEAN_NORMAL-temp$MEAN_CANCER) > fdr)
    wh <- intersect(wha, whb)
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    if(length(whx) != 0){

        tempy <- tempx[whx,]
        tempy$POSP <- tempy$POS/paired_sam[[2]][k]
        tempy$NEGP <- tempy$NEG/paired_sam[[2]][k]

        tpos <- c(tpos, tempy$POSP)
        tneg <- c(tneg, tempy$NEGP)
        tcancer <- c(tcancer, rep(all_cancer[k],length(tempy[[1]])))
        tdiff <- c(tdiff, tempy$MEAN_DIFF) ## MEAN_DIFF is mean ofthe paired normal and cancer differences
        tempn <- paste0(tempy$symbol,'_',tempy$as_id,'_',tempy$splice_type)
        asid <- c(asid, tempn)

    }

}

pdata <- data.frame(CANCER=tcancer, High=tpos, Low=tneg, impact=tdiff, ASID=asid)

##--- choose the max of the # of patients in which an event is gained or lost ----
flag <- c()
maxp <- c()
for(k in 1:length(pdata[[1]])){
    temp <- pdata[k,]
    maxp <- c(maxp, max(temp$High, temp$Low))
    if(temp$High > temp$Low){
        flag <- c(flag, 'High')
    }else{
        flag <- c(flag, 'Low')
    }
}

pdata$MAXP <- maxp
pdata$FLAG <- flag


## compute mean of max # of patients across as events in a cancer type
mdns <- c()
mxp <- c()
all_cancerx <- all_cancer[-4]
for(k in 1:length(all_cancerx)){
    tempx <- pdata[pdata$CANCER == all_cancerx[k], ]
    mdns <- c(mdns, signif(mean(tempx$MAXP)*100,3))
    mxp <- c(mxp, max(tempx$impact))
}
mdns <- paste0(mdns, '%')
pmdn <- data.frame(CANCER=all_cancerx, count=mdns,pos=mxp+0.25)

pdatax <- data.frame(cancer=tcancer, median_diff=tdiff)
p <- ggplot(pdatax, aes(cancer, median_diff)) + 
geom_jitter(aes(color=cancer),size=0.3)+
geom_violin(outlier.shape = NA)+
theme(legend.text=element_text(size=12))
basesize <- 8
p <- p + theme_classic(base_size = basesize) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="Mean \u0394PSI", limits=c(-1,1.2)) +
geom_text(data=pmdn, aes(y=pos, x=CANCER,label=count), position=position_dodge(width=0.9),hjust=0.5, vjust=0.5, angle=60, size=2.4)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#fb9a99','#ffff99',
    '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize), axis.title=element_text(size=basesize*1))+
guides(color='none')
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_dist.png"),width=3, height=2.8, dpi=500)


##-------- Pairwise overlap of slicing events affecting TFs across cancer types ----------------------------------------
all_files <- all_files[-4]
all_cancer <- substr(basename(all_files), 1,4)

TF_splicing_events <- list()
splicing_tfs <- list()
background_as <- c()

for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    wha <- which(temp$FDR < fdr)
    whb <- which(abs(temp$MEAN_NORMAL-temp$MEAN_CANCER) > fdr)
    wh <- intersect(wha, whb)
    tempx <- temp[wh, ]

    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    TF_splicing_events[[k]] <- tempy$as_id
    splicing_tfs[[k]] <- unique(tempy$symbol)
    background_as <- union(background_as, tempy$as_id)   
}


##--- Splicing events occuring in multiple cancer types --------------------------------
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
    temp_unn <- c()
    for(i in 1:loop1){
        temp_ovl <- TF_splicing_events[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- intersect(temp_ovl, TF_splicing_events[[which(all_cancer == temp_comb[[i]][j])]])
            }
        }
        temp_unn <- union(temp_unn, temp_ovl)
    }

    num_events_combs[[k]] <- temp_unn
    num_events_combs_counts <- c(num_events_combs_counts, length(temp_unn))
    cat('Cancer',k,'of',length(all_cancer),'done\n')
}

pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)

p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 8
maxv <- max(pdata$count)
p <- p + theme_classic(base_size = basesize) +
scale_x_discrete(name="Least # of cancer types") + 
scale_y_continuous(name="# of PTSEs", limits=c(0,maxv+300)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"),
panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
strip.text = element_text(size = basesize), axis.title=element_text(size=basesize))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Events_perturbing_overlap.png"),width=3, height=3, dpi=500)


##--- save excel file -----
tunn <- list()
tunnx <- list()

for(k in 1:length(num_events_combs[[1]])){
    temp_unn <- c()
    for(j in 1:length(all_cancer)){
        if(num_events_combs[[1]][k] %in% TF_splicing_events[[j]]){
            temp_unn <- c(temp_unn, all_cancer[j])
        }
    }
    tunn[[k]] <- paste(temp_unn, collapse=';')
    tunnx[[k]] <- temp_unn

}

wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Supplementary_Table_S2.xlsx'))
##--- column name explanations --
tdatat <- data.frame(`Table fields explanation`=c(
   'AS_ID: Alternative splicing ID provided by TCGA SpliceSeq',
     'CANCER_TYPES: Types of cancer (denoted as TCGA cancer codes) in which the considered splicing event in significantly perturbed'
     ))
openxlsx::addWorksheet(wb1, sheetName = 'INDEX')
openxlsx::writeData(wb1, sheet = 'INDEX', tdatat)
openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S2.xlsx'), overwrite = T)

tdatax <- data.frame(AS_ID=num_events_combs[[1]], CANCER_TYPES=unlist(tunn), len=lengths(tunnx))
tdatay <- tdatax[order(-tdatax$len), ]
tdatay <- tdatay[,-3]
openxlsx::addWorksheet(wb1, sheetName = 'PTSE-cancer combinations')
openxlsx::writeData(wb1, sheet = 'PTSE-cancer combinations', tdatay)
openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S2.xlsx'), overwrite = T)


##----- PSI profile of the pancancer (at least 10 cancer types) events -------------
events_to_consider <- num_events_combs[[10]]
tcancer <- c()
tvalue <- c()
tevent <- c()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    temp1 <- temp[temp$as_id %in% events_to_consider, ]
    temp1$ID <- paste0(temp1$symbol,'_',temp1$as_id,'_',temp1$splice_type)
    tcancer <- c(tcancer, rep(all_cancer[k], length(temp1[[1]])))
    tvalue <- c(tvalue, temp1$MEAN_DIFF)
    tevent <- c(tevent, temp1$ID)
}
pdata <- data.frame(tcancer=tcancer, val=tvalue, event=tevent)

pdata_mat <- as.data.frame(matrix(nrow=length(unique(pdata$event)), ncol=length(all_cancer),0))
rownames(pdata_mat) <- gtools::mixedsort(unique(pdata$event))
colnames(pdata_mat) <- gtools::mixedsort(unique(all_cancer))

for(k in 1:length(pdata[[1]])){
    wh1 <- which(rownames(pdata_mat) == pdata$event[k])
    wh2 <- which(colnames(pdata_mat) == pdata$tcancer[k])
    pdata_mat[wh1, wh2] <- as.numeric(pdata$val[k])
}

# pdata_mat <- cbind(data.frame(Event=unique(pdata$event)), pdata_mat)
pdx <- t(as.matrix(pdata_mat))#as.matrix(pdata_mat)#

p <- pheatmap(pdx,fontsize=3, cluster_rows=FALSE, cluster_cols=FALSE,cellheight=5, cellwidth = 5)
ggsave(p,filename=paste0(save_dir, "/Pancancer_events.png"),width=3.4, height=2, dpi=600)

saveRDS(num_events_combs, file = "data/overlapping_events.rds")

##-------------------------------------------------------------------------------


##----- unique splicing events -----------------------------------------------------
num_events_unq <- list()
num_tfs_unq <- list()

for(k in 1:length(all_cancer)){
    temp_unq <- TF_splicing_events[[k]]
    temp_un <- c()
    for(j in 1:length(all_cancer)){
        if(j != k){
            temp_un <- union(temp_un, TF_splicing_events[[j]])
        }
    }
    num_events_unq[[k]] <- setdiff(temp_unq, temp_un)

    temp_unq <- splicing_tfs[[k]]
    temp_un <- c()
    for(j in 1:length(all_cancer)){
        if(j != k){
            temp_un <- union(temp_un, splicing_tfs[[j]])
        }
    }
    num_tfs_unq[[k]] <- setdiff(temp_unq, temp_un)
}

pdata <- data.frame(cancer=all_cancer, count=lengths(num_events_unq))

saveRDS(num_events_unq, file = "data/unique_events.rds")

##--- Overlap with the survival associated events ---------------
## Study: OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers, NAR, 2022
## Downloaded http://47.98.127.64:8080/beta/download?fileName=spliceseq_clinical_as_survival.csv.gz
splicing <- data.table::fread('data/spliceseq_clinical_as_survival.csv')
splicing <- as.data.frame(splicing)

all_survival <- list()
for(k in 1:length(all_cancer)){
    temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
    all_survival[[k]] <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
}


##--- unique survival ---
survival_unq <- list()

for(k in 1:length(all_cancer)){
    temp_unq <- all_survival[[k]]
    temp_un <- c()
    for(j in 1:length(all_cancer)){
        if(j != k){
            temp_un <- union(temp_un, all_survival[[j]])
        }
    }
    survival_unq[[k]] <- setdiff(temp_unq, temp_un)
}


unq_events_surv <- list()

for(k in 1:length(all_cancer)){

    unq_events_surv[[k]] <- intersect(all_survival[[k]], num_events_unq[[k]])
    
}


pdata$Yes <- lengths(unq_events_surv)
pdata$No <- pdata$count-pdata$Yes
pdata1 <- pdata[,-2]

pdatax <- reshape2::melt(pdata1)
pdatax$vax <- rep(pdata$count, 2)
wh <- which(pdatax$variable == 'No')
pdatax$vax[wh] <- NA
p <- ggplot(pdatax, aes(cancer, value, fill=variable)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 8
p <- p + theme_classic(base_size = basesize) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of PTSEs", limits=c(0,max(pdata$count))) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
geom_text(aes(x=cancer, y=vax, label=value), position=position_stack(vjust=1), size=3)+
theme(axis.text.x = element_text(size = basesize, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"),
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),  
strip.text = element_text(size = basesize), axis.title=element_text(size=basesize), legend.position=c(0.8,0.86))+
guides(fill=guide_legend(title="Survival associated",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Sig_events_unique.png"),width=3.5, height=3, dpi=500)


##--- save excel file ---
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Supplementary_Table_S3.xlsx'))
##--- column name explanations --
tdatat <- data.frame(`Table fields explanation`=c(
   'AS_ID: Alternative splicing ID provided by TCGA SpliceSeq',
     'CI_Type: Type of clinical indicator',
     'Hazard_Ratio: Hazard ratio with respect to the considered clinical indicator'
     ))
openxlsx::addWorksheet(wb1, sheetName = 'INDEX')
openxlsx::writeData(wb1, sheet = 'INDEX', tdatat)
openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S3.xlsx'), overwrite = T)


temp_uq <- list()
for(k in 1:length(all_cancer)){

    temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
    all_surv <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
    surv <- rep(NA, length(num_events_unq[[k]]))
    citype <- rep(NA, length(num_events_unq[[k]]))

    for(j in 1:length(surv)){
        whx <- which(all_surv == num_events_unq[[k]][j])
        if(length(whx) > 0){
            surv[j] <- paste(temp1$Hazard_Ratio[whx], collapse=';')
            citype[j] <- paste(temp1$CI_Type[whx], collapse=';')
        }
    }

    tdatax <- data.frame(AS_ID=num_events_unq[[k]], CI_Type=citype, Hazard_Ratio=surv)
    temp_uq[[k]] <- tdatax
    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tdatax)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S3.xlsx'), overwrite = T)

}


##-- plot of unique splicing event hazard ratio and delta PSI values -----------
perturbed_tfs <- 'results_rep/TF_splicing/Supplementary_Table_S1.xlsx'
pdata <- data.frame(matrix(ncol=0,nrow=0))

for(k in 1:length(all_cancer)){
    kk <- k+1
    temp <- temp_uq[[k]]
    temp1_c <- temp[!is.na(temp$Hazard_Ratio), ]

    wh <- which(temp1_c[[2]] %like% ';')
    if(length(wh) != 0){
        temp1_c1 <- temp1_c[wh,]
        temp1_c2 <- temp1_c[-wh,]
        for(i in 1:length(temp1_c1[[1]])){
            tempv <- unlist(strsplit(temp1_c1[[2]][i], '[;]'))
            tempvv <- unlist(strsplit(temp1_c1[[3]][i], '[;]'))
            for(j in 1:length(tempv)){
                temp1_c2 <- rbind(temp1_c2, data.frame(AS_ID=temp1_c1$AS_ID[i], CI_Type=tempv[j], Hazard_Ratio=tempvv[j]))
            }
        }
    }else{
        temp1_c2 <- temp1_c
    }

    if(k < 4){
        tpsi <- openxlsx::read.xlsx(perturbed_tfs, kk)
    }else{
        tpsi <- openxlsx::read.xlsx(perturbed_tfs, kk+1)
    }


    pval <- c()
    gname <- c()
    ast <- c()
    as <- c()
    for(j in 1:length(temp1_c2[[1]])){
        pval <- c(pval, tpsi[tpsi$AS_ID == temp1_c2$AS_ID[j], ]$MEAN_DIFF)
        gname <- c(gname, tpsi[tpsi$AS_ID == temp1_c2$AS_ID[j], ]$SYMBOL)
        ast <- c(ast, tpsi[tpsi$AS_ID == temp1_c2$AS_ID[j], ]$AS_TYPE)
        as <- c(as, temp1_c2$AS_ID[j])
    }

    idd <- paste0(gname, '_', as, '_', ast)
    temp1_c2$MEAN_DIFF <- pval
    temp1_c2$CANCER <- rep(all_cancer[k], length(temp1_c2[[1]]))
    temp1_c2$ID <- idd
    pdata <- rbind(pdata, temp1_c2)
    
}

# whx <- which(abs(pdata$MEAN_DIFF) == max(abs(pdata$MEAN_DIFF)))
# whx <- which(abs(pdata$MEAN_DIFF) > 0.3)
pdata$IDD <- paste0(pdata$ID, '\n(',pdata$CI_Type, ')')
pdata$PL <- ""
pdata$PL[c(2, 117, 120, 200)] <- pdata$IDD[c(2, 117, 120, 200)]

pdata$Hazard_Ratio <- as.numeric(pdata$Hazard_Ratio)
pdata$MEAN_DIFF <- as.numeric(pdata$MEAN_DIFF)

pdata1 <- pdata[pdata$Hazard_Ratio < 3, ] ## filter outliers

p <- ggplot(pdata1, aes(Hazard_Ratio, MEAN_DIFF, color=CANCER, label=PL)) + 
# p <- ggplot(pdata, aes(AA, MEAN_DIFF, color=CANCER)) + 
geom_point()+
theme(legend.text=element_text(size=12))
basesize <- 8
p <- p + 
scale_x_continuous(name="Hazard ratio") + 
scale_y_continuous(name="Mean \u0394PSI", limits=c(-0.8, 0.8)) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#fb9a99','#ffff99',
    '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
geom_text_repel(family = "Poppins",
    max.overlaps=Inf,
                      size = 2.4,
                      color='black',
                      arrow = arrow(length = unit(0.010, "npc")),
                      min.segment.length = 0) +
theme_classic()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),axis.title=element_text(size=basesize),
    legend.text=element_text(size=basesize), legend.title=element_text(size=basesize),
    panel.border = element_blank())+
# guides(color='none')
guides(color=guide_legend(title="Cancer type", ncol=2, override.aes = list(size = 2)))
ggsave(p,filename=paste0(save_dir,"/Unique_HR.png"),width=5.5, height=2.8, dpi=500)


