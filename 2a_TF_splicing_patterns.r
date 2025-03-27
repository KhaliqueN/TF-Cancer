##############################################################################################
# Purpose: Plots to see 
# - The numebr of splicing events significantly different between normal and cancer tissues
# - Types of AS events in different cancer types
# - Number of AS events affecting TFs 
# - Types of AS events affecting TFs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

save_dir <- '../results'

## TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
##----------------------------------------------------------------------

input_dir <- '../data/PSI_data'
fdr <- 0.05
diff <- 0.1
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))

all_cancer <- substr(basename(all_files), 1,4)

num_sig_events_pos <- c()
num_sig_events_neg <- c()
samples <- c()
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()
RI <- c()

AAn <- c()
ADn <- c()
APn <- c()
ATn <- c()
ESn <- c()
MEn <- c()
RIn <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    ## save number of samples ----
    samples <- c(samples, length(which(colnames(temp) %like% 'TCGA')))
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    tempx_pos <- tempx[tempx$MEDIAN_DIFF > 0, ]
    tempx_neg <- tempx[tempx$MEDIAN_DIFF < 0, ]

    num_sig_events_pos <- c(num_sig_events_pos, length(tempx_pos[[1]]))
    num_sig_events_neg <- c(num_sig_events_neg, length(tempx_neg[[1]]))
    
    temp_count <- plyr::count(tempx_pos$splice_type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RI <- c(RI, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))

    temp_count <- plyr::count(tempx_neg$splice_type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AAn <- c(AAn, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    ADn <- c(ADn, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    APn <- c(APn, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    ATn <- c(ATn, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ESn <- c(ESn, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    MEn <- c(MEn, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RIn <- c(RIn, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))

}

##--- plot the number of splicing events -----------------
pdata1 <- data.frame(cancer=all_cancer, count=num_sig_events_pos, flag=rep('High in cancer'))
pdata2 <- data.frame(cancer=all_cancer, count=num_sig_events_neg, flag=rep('Low in cancer'))
pdata <- rbind(pdata1,pdata2)
p <- ggplot(pdata, aes(cancer, count, fill=flag)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(pdata1$count+pdata2$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of significant splicing events", limits=c(0,maxv+10)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')#
ggsave(p,filename=paste0(save_dir,"/Sig_events.png"),width=5.5, height=3, dpi=400)

##--- plot the number of splicing events of different types -----------------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=num_sig_events_pos, Cancer=all_cancer)
pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of splicing events")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_pos.png"),width=7, height=3, dpi=400)


pData <- data.frame(A=AAn, B=ADn, C=APn, D=ATn, E=ESn, F=MEn, G=RIn, X=num_sig_events_neg, Cancer=all_cancer)
pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of splicing events")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_neg.png"),width=7, height=3, dpi=400)



##---- Correlation with the numbers of perturbed events detected in existing studies
## Study: Alternative splicing perturbation landscape identifies RNA binding proteins as potential therapeutic targets in cancer, Molecular Therapy - Nucleic Acids, 2021
## Downloaded document 2, where we will take data from Supplementray Table S3
alien_file <- as.data.frame(readxl::read_excel('../data/1-s2.0-S2162253121000998-mmc2.xlsx', 3))
wh1 <- which(alien_file[[1]] == 'Up-regulated')
wh2 <- which(alien_file[[1]] == 'Down-regulated')
end1 <- wh2-1
strt2 <- wh2+1
upreg <- alien_file[wh1:end1, ]
upreg <- upreg[3:length(upreg[[1]]), ]
upregx <- rowSums(sapply(upreg[,-1], as.numeric ))
downreg <- alien_file[strt2:length(alien_file[[1]]), ]
downreg <- downreg[2:length(downreg[[1]]),]
downregx <- rowSums(sapply(downreg[,-1], as.numeric ))
alien_pert <- colSums(rbind(upregx, downregx))
wh <- which(upreg[[1]] %in% all_cancer)
pdatax <- data.frame(Cancer=all_cancer, mycount=num_sig_events_pos+num_sig_events_neg, aliencount=alien_pert[wh])

coral <- cor.test(x=pdatax$mycount, y=pdatax$aliencount, method = 'spearman')
p <- ggplot(pdatax, aes(mycount, aliencount, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 10
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of perturbed splicing events\n in this study", limits=c(0,max(pdatax$mycount))) + 
geom_text(aes(x=3000,y=11000, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,show.legend = FALSE)+
scale_y_continuous(name="# of perturbed splicing events\n in the Li et al. study", limits=c(0,max(pdatax$aliencount))) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0(save_dir,"/Sig_events_existing_Li.png"),width=4.5, height=2.5, dpi=400)


##---- Correlation with the numbers of perturbed events detected in existing studies
## Study: OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers, NAR, 2022
## Downloaded Supplementray Table S1 --> column 'SpliceSeq DASEs'
splicing <- readxl::read_excel('../data/Supplementary table 1.xlsx', 1)
splicing <- as.data.frame(splicing)
colnames(splicing) <- as.vector(splicing[1,])
splicing <- splicing[-1,]
splicing <- splicing[1:31,]
splicing1 <- splicing[,c(1,5,8,10)]
splicing1 <- splicing1[splicing1$`TCGA cancer` %in% all_cancer, ]
splicing1 <- splicing1[order(splicing1$`TCGA cancer`), ]
pdatax$aliencountZ <- as.numeric(splicing1[[4]])
wh <- which(!is.na(pdatax$aliencountZ))
if(length(wh) > 0){pdatax <- pdatax[wh,]}

coral <- cor.test(x=pdatax$mycount, y=pdatax$aliencountZ, method = 'spearman')
p <- ggplot(pdatax, aes(mycount, aliencountZ, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 10
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of perturbed splicing events\n in this study", limits=c(0,max(pdatax$mycount))) + 
geom_text(aes(x=2600,y=6000, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,show.legend = FALSE)+
scale_y_continuous(name="# of perturbed splicing events\n in the Zhang et al. study", limits=c(0,max(pdatax$aliencountZ))) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0(save_dir,"/Sig_events_existing_Zhang.png"),width=4.5, height=2.5, dpi=400)


##--- splicing events affecting TFs ----------------------------------------------------------------------
# num_tf_events <- c()
events_tf_pos <- list()
events_tf_neg <- list()
num_sig_events_pos <- c()
num_sig_events_neg <- c()
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()
RI <- c()

AAn <- c()
ADn <- c()
APn <- c()
ATn <- c()
ESn <- c()
MEn <- c()
RIn <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy_pos <- tempy[tempy$MEDIAN_DIFF > 0, ]
    tempy_neg <- tempy[tempy$MEDIAN_DIFF < 0, ]

    # num_tf_events <- c(num_tf_events, length(whx))
    events_tf_pos[[k]] <- tempy_pos$as_id
    events_tf_neg[[k]] <- tempy_neg$as_id

    num_sig_events_pos <- c(num_sig_events_pos, length(tempy_pos[[1]]))
    num_sig_events_neg <- c(num_sig_events_neg, length(tempy_neg[[1]]))
    
    temp_count <- plyr::count(tempy_pos$splice_type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RI <- c(RI, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))

    temp_count <- plyr::count(tempy_neg$splice_type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AAn <- c(AAn, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    ADn <- c(ADn, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    APn <- c(APn, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    ATn <- c(ATn, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ESn <- c(ESn, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    MEn <- c(MEn, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RIn <- c(RIn, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))
    
}

##--- plot the number of splicing events affecting TFs -----------------
pdata1 <- data.frame(cancer=all_cancer, count=num_sig_events_pos, flag=rep('High in cancer'))
pdata2 <- data.frame(cancer=all_cancer, count=num_sig_events_neg, flag=rep('Low in cancer'))
pdata <- rbind(pdata1,pdata2)
p <- ggplot(pdata, aes(cancer, count, fill=flag)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(pdata1$count+pdata2$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of significant splicing events", limits=c(0,maxv+10)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')#
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs.png"),width=5.5, height=3, dpi=400)


##--- plot the number of splicing events of different types affecting TFs -----------------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=num_sig_events_pos, Cancer=all_cancer)

pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of significant AS events \naffecting transcription factors")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_TFs_pos.png"),width=7, height=3, dpi=400)


pData <- data.frame(A=AAn, B=ADn, C=APn, D=ATn, E=ESn, F=MEn, G=RIn, X=num_sig_events_neg, Cancer=all_cancer)

pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of significant AS events \naffecting transcription factors")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_TFs_neg.png"),width=7, height=3, dpi=400)


##-------- Pairwise overlap of slicing events affecting TFs across cancer types --------------------------------------------
TF_splicing_events_pos <- list()
TF_splicing_events_neg <- list()
background_as_pos <- c()
background_as_neg <- c() 

for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
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
ggsave(p,filename=paste0(save_dir,'/Sig_events_TFs_overlap_pos.png'),width=4, height=2.5, dpi=400)


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
ggsave(p,filename=paste0(save_dir,'/Sig_events_TFs_overlap_neg.png'),width=4, height=2.5, dpi=400)


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
ggsave(p,filename=paste0(save_dir,'/Sig_events_TFs_overlap_pos_neg.png'),width=4, height=2.5, dpi=400)



##--- Splicing events occuring in multiple cancer types --------------------------------
combs <- list() ## store all combinations
for(k in 1:length(all_cancer)){
    combs[[k]] <- combn(all_cancer, k)
}

##-- speed up --read all files in memory!!
all_files_store <- list()
for(k in 1:length(all_files)){
    temp <- as.data.frame(data.table::fread(all_files[k], sep='\t'))
    ll <- length(temp) - 3
    all_files_store[[k]] <- read.table(all_files[k], colClasses = c('NULL','character','character',rep('NULL', ll)), header = TRUE, fill=TRUE)
}

num_events_combs_pos <- list()
num_events_combs_counts <- c()
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()
RI <- c()

for(k in 1:length(all_cancer)){

    temp_comb <- as.data.frame(combs[[k]])
    loop1 <- length(temp_comb)
    loop2 <- length(temp_comb[[1]])
    temp_unn <- c()
    tempy <- data.frame(matrix(ncol=0, nrow=0))
    for(i in 1:loop1){
        temp_ovl <- TF_splicing_events_pos[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- intersect(temp_ovl, TF_splicing_events_pos[[which(all_cancer == temp_comb[[i]][j])]])
            }
            temp <- as.data.frame(all_files_store[[which(all_cancer == temp_comb[[i]][2])]])
        }
        scols <- which(colnames(temp) %in% c('as_id','splice_type'))
        tempx <- temp[temp$as_id %in% temp_ovl, scols]
        tempy <- rbind(tempy, tempx)
        temp_unn <- union(temp_unn, temp_ovl)
    }

    tempy <- unique(tempy)
    num_events_combs_pos[[k]] <- temp_unn
    num_events_combs_counts <- c(num_events_combs_counts, length(temp_unn))

    temp_count <- plyr::count(tempy$splice_type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RI <- c(RI, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))
    cat('Cancer',k,'of',length(all_cancer),'done\n')

}

pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of significant splicing events", limits=c(0,(max(pdata$count))+500)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_overlap_pos.png"),width=7, height=3, dpi=400)


##---
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=num_events_combs_counts, Cancer=as.factor(seq(1,length(all_cancer))))

pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("# of cancer types")+ylab("% of splicing events")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_TFs_overlap_pos.png"),width=7, height=3, dpi=400)



num_events_combs_neg <- list()
num_events_combs_counts <- c()
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()
RI <- c()

for(k in 1:length(all_cancer)){

    temp_comb <- as.data.frame(combs[[k]])
    loop1 <- length(temp_comb)
    loop2 <- length(temp_comb[[1]])
    temp_unn <- c()
    tempy <- data.frame(matrix(ncol=0, nrow=0))
    for(i in 1:loop1){
        temp_ovl <- TF_splicing_events_neg[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- intersect(temp_ovl, TF_splicing_events_neg[[which(all_cancer == temp_comb[[i]][j])]])
            }
            temp <- as.data.frame(all_files_store[[which(all_cancer == temp_comb[[i]][2])]])
        }
        scols <- which(colnames(temp) %in% c('as_id','splice_type'))
        tempx <- temp[temp$as_id %in% temp_ovl, scols]
        tempy <- rbind(tempy, tempx)
        temp_unn <- union(temp_unn, temp_ovl)
    }

    tempy <- unique(tempy)
    num_events_combs_neg[[k]] <- temp_unn
    num_events_combs_counts <- c(num_events_combs_counts, length(temp_unn))

    temp_count <- plyr::count(tempy$splice_type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RI <- c(RI, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))
    cat('Cancer',k,'of',length(all_cancer),'done\n')

}

pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of significant splicing events", limits=c(0,(max(pdata$count))+500)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_overlap_neg.png"),width=7, height=3, dpi=400)


##---
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=num_events_combs_counts, Cancer=as.factor(seq(1,length(all_cancer))))

pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("# of cancer types")+ylab("% of splicing events")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_TFs_overlap_neg.png"),width=7, height=3, dpi=400)



##----- PSI profile of the pancancer (at least 10 cancer types) events -------------
events_to_consider <- num_events_combs_pos[[10]]
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

pdata1 <- data.frame(tcancer=tcancer, val=tvalue, event=tevent)

events_to_consider <- num_events_combs_neg[[10]]
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

p <- pheatmap(pdx)
ggsave(p,filename=paste0(save_dir, "/Pancancer_events.png"),width=7, height=5, dpi=400)


##----- unique splicing events -----------------------------------------------------
num_events_unq <- list()
num_events_unq_counts <- c()
unique_events_tf <- list()
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()
RI <- c()

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

    temp <- data.table::fread(all_files[k], sep='\t')
    tempy <- temp[which(temp$as_id %in% temp_unq), ]
    unique_events_tf[[k]] <- tempy$as_id
    temp_count <- plyr::count(tempy$splice_type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RI <- c(RI, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))

}

pdata <- data.frame(cancer=all_cancer, count=num_events_unq_counts)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of significant splicing events", limits=c(0,150)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Unique_sig_events_TFs_pos.png"),width=3.5, height=3, dpi=400)

##------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=num_events_unq_counts, Cancer=all_cancer)

pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of splicing events \naffecting transcription factors")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Unique_sig_events_types_TFs_pos.png"),width=7, height=3, dpi=400)



num_events_unq <- list()
num_events_unq_counts <- c()
unique_events_tf <- list()
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()
RI <- c()

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

    temp <- data.table::fread(all_files[k], sep='\t')
    tempy <- temp[which(temp$as_id %in% temp_unq), ]
    unique_events_tf[[k]] <- tempy$as_id
    temp_count <- plyr::count(tempy$splice_type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RI <- c(RI, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))

}

pdata <- data.frame(cancer=all_cancer, count=num_events_unq_counts)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of significant splicing events", limits=c(0,150)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Unique_sig_events_TFs_neg.png"),width=3.5, height=3, dpi=400)

##------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=num_events_unq_counts, Cancer=all_cancer)

pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of splicing events \naffecting transcription factors")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Unique_sig_events_types_TFs_neg.png"),width=7, height=3, dpi=400)



##--- work from here -----











##--- Correlation with the number of samples --------------------
pData$Samples <- samples
coral <- cor.test(x=pData$X, y=pData$Samples, method = 'spearman')
p <- ggplot(pData, aes(X, Samples, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 10
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="% of significant AS events \naffecting transcription factors", limits=c(0,max(pData$X))) + 
geom_text(aes(x=300,y=900, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,show.legend = FALSE)+
scale_y_continuous(name="# of samples", limits=c(0,max(pData$Samples))) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_samples.png"),width=4.5, height=2.5, dpi=400)


##--- Overlap with the survival associated events ---------------
## Study: OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers, NAR, 2022
## Downloaded http://47.98.127.64:8080/beta/download?fileName=spliceseq_clinical_as_survival.csv.gz
splicing <- data.table::fread('../data/spliceseq_clinical_as_survival.csv')
splicing <- as.data.frame(splicing)
events_tf_surv <- list()
events_tf_surv_c <- c()
events_tf_surv_p <- c()
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()
RI <- c()

OS <- c()
DFI <- c()
PFI <- c()
DSS <- c()

for(k in 1:length(all_cancer)){
    temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
    temp1$SE <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
    whi <- intersect(temp1$SE, events_tf[[k]])
    tempy <- temp1[temp1$SE %in% whi, ]

    # # temp1 <- temp[temp$CI_Type == 'overall_survival', ]$Splice_Event
    # temp2 <- unlist(lapply(strsplit(temp1, '[_]'), '[[', 3))
    events_tf_surv[[k]] <- intersect(tempy$SE, events_tf[[k]])
    events_tf_surv_c <- c(events_tf_surv_c, length(events_tf_surv[[k]]))
    events_tf_surv_p <- c(events_tf_surv_p, (length(events_tf_surv[[k]])/length(events_tf[[k]])*100))

    temp_count <- plyr::count(tempy$Splice_Type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RI <- c(RI, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))

    temp_count <- plyr::count(tempy$CI_Type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    OS <- c(OS, ifelse(length(which(temp_count$x == 'overall_survival') != 0),temp_count$prct[which(temp_count$x == 'overall_survival')],0))
    PFI <- c(PFI, ifelse(length(which(temp_count$x == 'progression_free_interval') != 0),temp_count$prct[which(temp_count$x == 'progression_free_interval')],0))
    DFI <- c(DFI, ifelse(length(which(temp_count$x == 'disease_free_interval') != 0),temp_count$prct[which(temp_count$x == 'disease_free_interval')],0))
    DSS <- c(DSS, ifelse(length(which(temp_count$x == 'disease_specific_survival') != 0),temp_count$prct[which(temp_count$x == 'disease_specific_survival')],0))
    
}

pdata <- data.frame(cancer=all_cancer, count=events_tf_surv_c, pcount=events_tf_surv_p)
p <- ggplot(pdata, aes(cancer, pcount)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="% of events affecting TFs\n that are survival associated", limits=c(0,75)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_surv.png"),width=3.5, height=3, dpi=400)

# ##--- plot the number of survival associated splicing events of different types affecting TFs -----------------
# pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=events_tf_surv_c, Cancer=all_cancer)
# pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
# pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
# basesize <- 8
# ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
# geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
# scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
# xlab("Cancer type")+ylab("% of significant AS events \naffecting transcription factors")+
# scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
#     "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
# values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
# theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
#     axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
#     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"))+
# guides(fill=guide_legend(title="Alternative splicing type"))
# ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_TFs_surv.png"),width=7, height=3, dpi=400)


# ##--- plot the number of survival associated splicing events of different survival types -----------------
# pData <- data.frame(A=OS, B=PFI, C=DFI, D=DSS, X=events_tf_surv_c, Cancer=all_cancer)
# pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
# pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
# basesize <- 8
# ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
# geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
# scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
# xlab("Cancer type")+ylab("% of significant AS events \naffecting transcription factors")+
# scale_fill_manual(labels=c("A" = "Overall survival", 
#     "B"="Progression free interval", "C"="Disease free interval","D"="Disease specific survival"), 
# values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
# theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
#     axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
#     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"))+
# guides(fill=guide_legend(title="Alternative splicing type"))
# ggsave(ppx,filename=paste0(save_dir, "/Surv_types_TFs_surv.png"),width=7, height=3, dpi=400)


gene_data <- list(TF_splicing_events, num_events_combs, num_events_unq)
saveRDS(gene_data, '../data/TF_splicing_events_cancer.rds')


##--- Unique splicing evevnts overlap with the survival associated events ---------------
## Study: OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers, NAR, 2022
## Downloaded http://47.98.127.64:8080/beta/download?fileName=spliceseq_clinical_as_survival.csv.gz
splicing <- data.table::fread('../data/spliceseq_clinical_as_survival.csv')
splicing <- as.data.frame(splicing)
events_tf_surv <- list()
events_tf_surv_c <- c()
events_tf_surv_p <- c()
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()
RI <- c()

OS <- c()
DFI <- c()
PFI <- c()
DSS <- c()

for(k in 1:length(all_cancer)){
    temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
    temp1$SE <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
    whi <- intersect(temp1$SE, unique_events_tf[[k]])
    tempy <- temp1[temp1$SE %in% whi, ]

    # # temp1 <- temp[temp$CI_Type == 'overall_survival', ]$Splice_Event
    # temp2 <- unlist(lapply(strsplit(temp1, '[_]'), '[[', 3))
    events_tf_surv[[k]] <- intersect(tempy$SE, unique_events_tf[[k]])
    events_tf_surv_c <- c(events_tf_surv_c, length(events_tf_surv[[k]]))
    events_tf_surv_p <- c(events_tf_surv_p, (length(events_tf_surv[[k]])/length(unique_events_tf[[k]])*100))

    temp_count <- plyr::count(tempy$Splice_Type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    RI <- c(RI, ifelse(length(which(temp_count$x == 'RI') != 0),temp_count$prct[which(temp_count$x == 'RI')],0))

    temp_count <- plyr::count(tempy$CI_Type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    OS <- c(OS, ifelse(length(which(temp_count$x == 'overall_survival') != 0),temp_count$prct[which(temp_count$x == 'overall_survival')],0))
    PFI <- c(PFI, ifelse(length(which(temp_count$x == 'progression_free_interval') != 0),temp_count$prct[which(temp_count$x == 'progression_free_interval')],0))
    DFI <- c(DFI, ifelse(length(which(temp_count$x == 'disease_free_interval') != 0),temp_count$prct[which(temp_count$x == 'disease_free_interval')],0))
    DSS <- c(DSS, ifelse(length(which(temp_count$x == 'disease_specific_survival') != 0),temp_count$prct[which(temp_count$x == 'disease_specific_survival')],0))
    
}

pdata <- data.frame(cancer=all_cancer, count=events_tf_surv_c, pcount=events_tf_surv_p)
p <- ggplot(pdata, aes(cancer, pcount)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="% of events affecting TFs\n that are survival associated", limits=c(0,40)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Unique_sig_events_TFs_surv.png"),width=3.5, height=3, dpi=400)

##--- plot the number of survival associated splicing events of different types affecting TFs -----------------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=events_tf_surv_c, Cancer=all_cancer)
pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of significant AS events \naffecting transcription factors")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Unique_sig_events_types_TFs_surv.png"),width=7, height=3, dpi=400)


##--- plot the number of survival associated splicing events of different survival types -----------------
pData <- data.frame(A=OS, B=PFI, C=DFI, D=DSS, X=events_tf_surv_c, Cancer=all_cancer)
pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","G","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of significant AS events \naffecting transcription factors")+
scale_fill_manual(labels=c("A" = "Overall survival", 
    "B"="Progression free interval", "C"="Disease free interval","D"="Disease specific survival"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Unique_surv_types_TFs_surv.png"),width=7, height=3, dpi=400)



