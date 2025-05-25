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
num_sams <- 0.5
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

wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Perturbed_TF_splicing_events.xlsx'))

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]

    ## save excel sheet ----
    tempz <- data.frame(tempy)
    tdatat <- tempz[, c(1,2,3, seq(length(tempz)-9, length(tempz)))]
    tdatat <- tdatat[order(-abs(tdatat$MEDIAN_DIFF)), ]

    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tdatat)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Perturbed_TF_splicing_events.xlsx'), overwrite = T)

    events_tf[[k]] <- tempy$as_id
    
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

##--- plot the number of splicing events affecting TFs -----------------
pdata <- data.frame(cancer=all_cancer, count=lengths(events_tf))
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(pdata$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of significant splicing events", limits=c(0,maxv+120)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs.png"),width=3.5, height=3, dpi=400)

##--- plot the number of splicing events of different types affecting TFs -----------------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, G=RI, X=lengths(events_tf), Cancer=all_cancer)
pData <- pData[pData$Cancer != 'ESCA', ]
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
ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_TFs.png"),width=7, height=3, dpi=400)



###------ Distribution of perturbation values ----------------------------------------------------
##------------------------------------------------
tpos <- c()
tneg <- c()
tcancer <- c()
tdiff <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy$POSP <- tempy$POS/paired_sam[[2]][k]
    tempy$NEGP <- tempy$NEG/paired_sam[[2]][k]

    tpos <- c(tpos, tempy$POSP)
    tneg <- c(tneg, tempy$NEGP)
    tcancer <- c(tcancer, rep(all_cancer[k],length(tempy[[1]])))
    tdiff <- c(tdiff, tempy$MEDIAN_DIFF)

}

pdata <- data.frame(CANCER=tcancer, High=tpos, Low=tneg, impact=tdiff)
# p <- ggplot(pdata, aes(High, impact)) + 
# geom_point(size=0.2)+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_y_continuous(name="Difference between \nmedian PSIs of paired samples", limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) + 
# scale_x_continuous(name="Fraction of patients in a cancer type") +
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# geom_hline(yintercept=0.2, color='red', linetype='dashed')+
# geom_hline(yintercept=-0.2, color='red', linetype='dashed')+
# guides(color='none')
# ggsave(p,filename=paste0(save_dir,"/Gained_events_patients.png"),width=3.5, height=3, dpi=400)

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

p <- ggplot(pdata, aes(MAXP, impact, color=FLAG)) + 
geom_point(size=0.2)+geom_density_2d(colour='white', size=0.2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_y_continuous(name="Difference between \nmedian PSIs of paired samples", limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) + 
scale_x_continuous(name="Fraction of patients in a cancer type", limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# geom_hline(yintercept=0.2, color='black', linetype='dashed', size=0.2)+
# geom_hline(yintercept=-0.2, color='black', linetype='dashed', size=0.2)+
scale_color_manual(values=c('#e41a1c','#377eb8'))+
guides(color=guide_legend(title="PSI in cancer\nvs.\npaired normal\nsamples",ncol=1, override.aes = list(size=2)))
ggsave(p,filename=paste0(save_dir,"/Perturbed_events_patients.png"),width=4.5, height=3, dpi=600)


pdata <- data.frame(cancer=tcancer, median_diff=tdiff)
p <- ggplot(pdata, aes(cancer, median_diff)) + 
geom_jitter(aes(color=cancer),size=0.2)+
geom_violin(outlier.shape = NA)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="Difference between \nmedian PSIs of paired samples", limits=c(-1,1)) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99',
    '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color='none')
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_dist.png"),width=5, height=3, dpi=400)

###---------------------------------------------------------------------------------------------------------------------





##-------- Pairwise overlap of slicing events affecting TFs across cancer types ----------------------------------------
all_files <- all_files[-4]
all_cancer <- substr(basename(all_files), 1,4)

TF_splicing_events <- list()
background_as <- c()

for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]

    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    TF_splicing_events[[k]] <- tempy$as_id
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









##--- Number of TFs with perturbed splicing events ----------------------------------------------------------------------
num_tf_events <- c()
tf_events <- list()
for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    num_tf_events <- c(num_tf_events, length(unique(tempy$symbol)))
    tf_events[[k]] <- unique(tempy$symbol)
    
}

##--- plot the number of splicing events affecting TFs -----------------
pdata <- data.frame(cancer=all_cancer, count=num_tf_events)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of TFs affected by \nperturbed splicing events", limits=c(0,max(pdata$count+50))) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_genes.png"),width=3.5, height=3, dpi=400)

