##############################################################################################
# Purpose: Plots to see 
# - The numebr of EVENTS significantly different between normal and cancer tissues
# - Types of AS events in different cancer types
# - Number of AS events affecting TFs 
# - Types of AS events affecting TFs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)
library(dplyr)

save_dir <- '../results_new/TF_splicing'

if(dir.exists(save_dir)){
	unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

## TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
paired_sam <- data.table::fread('../data/cancer_paired_samples.txt')
##----------------------------------------------------------------------

input_dir <- '../data/PSI_data'
fdr <- 0.05
# num_sams <- 0.5
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
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

wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Perturbed_TF_splicing_events.xlsx'))

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

    colnames(tdatat) <- c('SYMBOL', 'AS_ID', 'SPLICE_TYPE', 'P_VALUE', 'FDR',
    'MEDIAN_DIFF', 'MEDIAN_NORMAL', 'MEDIAN_CANCER', 'MEAN_DIFF', 'MEAN_NORMAL', 
    'MEAN_CANCER', 'CANCER_HIGH', 'NORMAL_HIGH', 'PAIRED_SAMPLES' )
    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tdatat)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Perturbed_TF_splicing_events.xlsx'), overwrite = T)

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
pdatae <- data.frame(cancer=all_cancer, count1=lengths(events_tf), count2=lengths(num_tfs))
maxv <- max(pdatae$count1)
colnames(pdatae) <- c('Cancer','PTSE','TF')
pdata <- reshape2::melt(pdatae)
pdata[pdata == 0] <- NA
p <- ggplot(pdata, aes(Cancer, value, fill=variable)) + 
geom_bar(stat="identity",position="dodge")+
theme(legend.text=element_text(size=12))
basesize <- 8
p <- p + theme_bw(base_size = basesize) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of PTSEs or TFs", limits=c(0,maxv+100)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
scale_fill_manual(values=c('#d95f02','#1b9e77'))+
geom_text(aes(label=value), position=position_dodge(width=0.9),hjust=0, vjust=0.5, angle=85, size=3)+
theme(axis.text.x = element_text(size = basesize, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize), axis.title=element_text(basesize*1.25), legend.position=c(0.85,0.82))+
guides(fill=guide_legend(title="Entity",ncol=1))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs.png"),width=3.5, height=3, dpi=400)


# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity")+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# maxv <- max(pdata$count)
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of significant splicing events", limits=c(0,maxv+120)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill=guide_legend(title="Percent Spliced In",ncol=1))
# ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs.png"),width=3.5, height=3, dpi=400)


# ##--- percentage of perturbed TF splicing events with respect to background ----
# bg_tfs <- lengths(num_all_events_tfs)/lengths(num_all_events)
# sig_tfs <- lengths(num_sig_events_tfs)/lengths(num_sig_events)

# ##------------------------------------------------------------------------------

# ##--- correlation with number of paired samples ------------------------
# pdatac <- data.frame(cancer=all_cancer, count=paired_sam$PAIRED_SAMPLES)
# pdatac$events <- pdatae$PTSEs
# pdatac <- pdatac[complete.cases(pdatac),]
# coral <- cor.test(x=pdatac$events, y=pdatac$count, method = 'spearman')

#     p <- ggplot(pdata, aes(Perturbed, Survival, color=Cancer)) + 
#     geom_point(size=2)+
#     theme(legend.text=element_text(size=12))
#     basesize <- 12
#     p <- p + theme_bw(base_size = basesize * 0.8) +
#     scale_x_continuous(name="# of perturbed edges") + 
#     scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
#     geom_text(aes(x=x3[ntype],y=y3[ntype], label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
#     scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
#     theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
#     axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
#     strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
#     guides(color=guide_legend(title="Cancer type",ncol=3))
#     ggsave(p,filename=paste0(save_dir,"/",net_type[ntype],"_correlation_perturbed.png"),width=4.8, height=2.5, dpi=400)


# ##--------------------------------------------------

##--- plot the number of splicing events of different types affecting TFs -----------------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=lengths(events_tf), Cancer=all_cancer)
pData <- pData[pData$Cancer != 'ESCA', ]
pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
# geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,101), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of PTSEs")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"), legend.position="bottom",axis.title=element_text(size=basesize),
    legend.text=element_text(size=basesize), legend.title=element_text(size=basesize))+
guides(fill=guide_legend(title="Type of\nalternative\nsplicing event", ncol=1, override.aes = list(size = 1)))
ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_TFs.png"),width=3, height=4, dpi=600)




###----- UMAP/TSNE ------------------------------------------------------------------------------------
casidx <- tumap[[1]]
for(k in 2:length(all_cancer)){
    casidx <- intersect(casidx, tumap[[k]])
}

alltfe <- events_tf[[1]]
for(k in 2:length(all_cancer)){
   alltfe <- union(alltfe, events_tf[[k]])
}

casid <- intersect(casidx, alltfe)


umap_data <- data.frame(matrix(ncol=0, nrow=0))
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
    umap_data <- rbind(umap_data, tempy)
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
basesize <- 10
p <- p + theme_bw(base_size = basesize) +
scale_x_continuous(name="tSNE dimension 1") + 
scale_y_continuous(name="tSNE dimension 2") +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99',
    '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# geom_text(aes(label=value), position=position_stack(vjust=0.5), size=3)+
theme(axis.text.x = element_text(size = basesize*1.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize*1.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize*1.6), axis.title=element_text(size=basesize*1.6), legend.position="bottom", 
legend.text=element_text(size=basesize*1.3), legend.title=element_text(size=basesize*1.3))+
guides(color=guide_legend(title="Cancer\ntype", ncol=4, override.aes = list(size = 3)))
ggsave(p,filename=paste0(save_dir,"/tSNE.png"),width=5, height=6, dpi=500)


# basesize <- 10
# p <- ggplot(umap_reduction_df, aes(V1, V2, color=Sample)) + 
# geom_point(alpha=0.75, size=0.2)+facet_wrap(~Cancer,ncol=5)+
# theme(legend.text=element_text(size= basesize * 0.8))
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_continuous(name="tSNE dimension 1") + 
# scale_y_continuous(name="tSNE dimension 2") +
# scale_color_manual(values=c('#e31a1c','#1f78b4'))+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(color=guide_legend(title="Sample type", ncol=1, override.aes = list(size = 2)))
# ggsave(p,filename=paste0(save_dir,"/tSNE_all.png"),width=6, height=3, dpi=600)


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
        tdiff <- c(tdiff, tempy$MEAN_DIFF)
        tempn <- paste0(tempy$symbol,'_',tempy$as_id,'_',tempy$splice_type)
        asid <- c(asid, tempn)

    }

}

pdata <- data.frame(CANCER=tcancer, High=tpos, Low=tneg, impact=tdiff, ASID=asid)

##--- choose the max of # of patients in which an event is gained or lost ----
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

p <- ggplot(pdata, aes(MAXP, impact, color=CANCER)) + 
geom_point(alpha=0.6)+geom_density_2d(colour='white', size=0.2)+
theme(legend.text=element_text(size=12))
basesize <- 8
p <- p + theme_bw(base_size = basesize) +
scale_y_continuous(name="Mean \u0394PSI (cancer - normal) \nof paired samples", limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) + 
scale_x_continuous(name="# of patients in a cancer type", limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
theme(axis.text.x = element_text(size = basesize, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize , angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize ), axis.title=element_text(basesize))+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#fb9a99','#ffff99',
    '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
guides(color=guide_legend(title="Cancer type",ncol=2, override.aes = list(size=2)))
ggsave(p,filename=paste0(save_dir,"/Perturbed_events_patients.png"),width=5, height=3.5, dpi=600)


# md_pdata <- pdata[,c(1,6)]
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
p <- p + theme_bw(base_size = basesize) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="Mean \u0394PSI", limits=c(-1,1.2)) +
geom_text(data=pmdn, aes(y=pos, x=CANCER,label=count), position=position_dodge(width=0.9),hjust=0.5, vjust=0.5, angle=60, size=2.5)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#fb9a99','#ffff99',
    '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize), axis.title=element_text(basesize*1.25))+
guides(color='none')
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_dist.png"),width=3.5, height=2.5, dpi=400)




##-------- Pairwise overlap of slicing events affecting TFs across cancer types ----------------------------------------
all_files <- all_files[-4]
all_cancer <- substr(basename(all_files), 1,4)

TF_splicing_events <- list()
splicing_tfs <- list()
background_as <- c()

for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
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
basesize <- 12
maxv <- max(pdata$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of splicing events", limits=c(0,maxv+300)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"),
panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Events_perturbing_overlap.png"),width=3.5, height=3, dpi=300)


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
    tvalue <- c(tvalue, temp1$MEDIAN_DIFF)
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
pdx <- t(as.matrix(pdata_mat))

p <- pheatmap(pdx,fontsize=3, cluster_rows=FALSE, cluster_cols=FALSE,cellheight=5, cellwidth = 5)
ggsave(p,filename=paste0(save_dir, "/Pancancer_events.png"),width=7, height=5, dpi=600)


# ##--- Number of TFs with perturbed splicing events ----------------------------------------------------------------------
# num_tf_events <- c()
# tf_events <- list()
# for(k in 1:length(all_cancer)){

#     temp <- data.table::fread(all_files[k], sep='\t')
#     wh <- which(temp$FDR < fdr)
#     tempx <- temp[wh, ]
#     whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
#     tempy <- tempx[whx,]
#     num_tf_events <- c(num_tf_events, length(unique(tempy$symbol)))
#     tf_events[[k]] <- unique(tempy$symbol)
    
# }

# ##--- plot the number of splicing events affecting TFs -----------------
# pdata <- data.frame(cancer=all_cancer, count=num_tf_events)
# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity",position=position_dodge())+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of TFs affected by \nperturbed splicing events", limits=c(0,max(pdata$count+50))) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
# ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_genes.png"),width=3.5, height=3, dpi=400)




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


##--- Overlap with the survival associated events ---------------
## Study: OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers, NAR, 2022
## Downloaded http://47.98.127.64:8080/beta/download?fileName=spliceseq_clinical_as_survival.csv.gz
splicing <- data.table::fread('../data/spliceseq_clinical_as_survival.csv')
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

p <- ggplot(pdatax, aes(cancer, value, fill=variable)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of events", limits=c(0,max(pdata$count))) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
geom_text(aes(label=value), position=position_stack(vjust=0.5), size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Survival",ncol=1))
ggsave(p,filename=paste0(save_dir,"/Sig_events_unique.png"),width=5, height=3, dpi=400)


# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity")+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# maxv <- max(pdata$count)
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of significant splicing events", limits=c(0,maxv+120)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill=guide_legend(title="Percent Spliced In",ncol=1))
# ggsave(p,filename=paste0(save_dir,"/Sig_events_unique.png"),width=3.5, height=3, dpi=400)


# pdata <- data.frame(cancer=all_cancer, count=events_tf_surv_c, pcount=events_tf_surv_p)
# p <- ggplot(pdata, aes(cancer, pcount)) + 
# geom_bar(stat="identity",position=position_dodge())+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="% of events affecting TFs\n that are survival associated", limits=c(0,75)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
# ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_surv.png"),width=3.5, height=3, dpi=400)




