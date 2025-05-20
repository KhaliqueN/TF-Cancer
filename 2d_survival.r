##############################################################################################
# Purpose: Survival associated events
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

save_dir <- '../results_new/Survival'

if(!dir.exists(save_dir)){
	dir.create(save_dir, recursive=TRUE)
}

##-- TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
##----------------------------------------------------------------------

input_dir <- '../data/PSI_data'
fdr <- 0.05
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))
all_cancer <- substr(basename(all_files), 1,4)


##------------------------------------------------
tf_events <- list()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[which(all_files %like% all_cancer[k])], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy <- tempy[order(-abs(tempy$MEDIAN_DIFF)), ]
    tf_events[[k]] <- unique(tempy$as_id) 
}



##--- Overlap with the survival associated events ---------------
## Study: OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers, NAR, 2022
## Downloaded http://47.98.127.64:8080/beta/download?fileName=spliceseq_clinical_as_survival.csv.gz
splicing <- data.table::fread('../data/spliceseq_clinical_as_survival.csv')
splicing <- as.data.frame(splicing)
events_tf_surv <- list()
events_tf_surv_sig <- list()

events_tf_surv_c <- c()
tf_hr <- c()
tf_medians <- c()
tcancer <- c()
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Perturbed_TF_splicing_events_survival.xlsx'))

for(k in 1:length(all_cancer)){

    events_tfx <- tf_events[[k]]
    temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
    temp1$SE <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
    whi <- intersect(temp1$SE, events_tfx)
    tempy <- temp1[temp1$SE %in% whi, ]

    ##--
    tempsi <- data.table::fread(all_files[k])
    tempsi <- tempsi[tempsi$as_id %in% whi, ]

    tmdif <- c()
    tfdr <- c()

    for(j in 1:length(tempy[[1]])){
      tmdif <- c(tmdif, tempsi$MEDIAN_DIFF[which(tempsi$as_id == tempy$SE[j])])
      tfdr <- c(tfdr, tempsi$FDR[which(tempsi$as_id == tempy$SE[j])])
    }

    tempy$MEDIAN_DIFF <- tmdif
    tempy$FDR <- tfdr
    tempy <- tempy[order(-abs(tempy$MEDIAN_DIFF)),]
    tempz <- tempy[tempy$MEDIAN_DIFF > 0.2,] 

    events_tf_surv[[k]] <- intersect(tempy$SE, events_tfx)
    events_tf_surv_sig[[k]] <- intersect(tempz$SE, events_tfx)

    events_tf_surv_c <- c(events_tf_surv_c, length(events_tf_surv[[k]]))

    ##--- add HR values -------------------------------------
    tf_hr <- c(tf_hr, tempy$Hazard_Ratio)
    tf_medians <- c(tf_medians, tempy$MEDIAN_DIFF)
    tcancer <- c(tcancer, rep(all_cancer[k], length(tempy$Hazard_Ratio)))
    ##-------------------------------------------------------

    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tempy)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Perturbed_TF_splicing_events_survival.xlsx'), overwrite = T)

}

pdata <- data.frame(cancer=all_cancer, count=events_tf_surv_c)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="% of survival associated \nperturbed events", limits=c(0,max(pdata$count)+50)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_surv.png"),width=3.5, height=3, dpi=400)


pdata <- data.frame(cancer=tcancer, HR=log2(tf_hr), MEDIANDIFF=tf_medians)
wh1 <- which(pdata$HR > 10)
wh2 <- which(pdata$HR < -10)
wh <- union(wh1, wh2)
pdata1 <- pdata[-wh,]
basesize <- 10
ppx <- ggplot(data = pdata1, aes(x=HR, y=MEDIANDIFF)) + 
geom_point(alpha=0.6)+
scale_x_continuous()+
scale_y_continuous()+
geom_hline(yintercept=0.2, color='red', linetype='dashed')+
geom_hline(yintercept=-0.2, color='red', linetype='dashed')+
xlab("Hazard ratio (log2)")+ylab("Difference between \nmedians of paired samples")+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill='none')
ggsave(ppx,filename=paste0(save_dir, "/Scatter_hazardRatio.png"),width=4, height=3, dpi=400)



##--- survival plots -----------------------------

for(k in 1:length(all_cancer)){

    temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
    temp1$SE <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
    temp <- temp1[which(temp1$SE %in% events_tf_surv_sig[[k]]),]

    if(nrow(temp) == 0){next}

    tempx <- temp[,c(1,14,17)]

    ids <- strsplit(tempx[[1]], '[_]')
    tgene <- unlist(lapply(ids, '[[', 1))
    tspt <- unlist(lapply(ids, '[[', 2))
    tevt <- unlist(lapply(ids, '[[', 3))

    temp_event <- data.table::fread(all_files[which(all_files %like% all_cancer[k])], sep='\t')
    temp_ids <- c()
    mval <- c()
    for(j in 1:length(tevt)){
        mval <- c(mval, temp_event[temp_event$as_id == tevt[j], ]$MEDIAN_DIFF)
        temp_ids <- c(temp_ids, paste0(tgene[j], '_', tevt[j],'_', tspt[j], ' (',temp_event[temp_event$as_id == tevt[j], ]$MEDIAN_DIFF,')'))
    }

    tempx$Splice_Event <- temp_ids
    tempx$HR <- log2(tempx$Hazard_Ratio)
    tempx$MEDIAN_DIFF <- mval
    tempx <- tempx[order(-tempx$MEDIAN_DIFF), ]
    tempx$CI_Type <- gsub('_','\n',tempx$CI_Type)
    tempx$Splice_Event <- factor(tempx$Splice_Event, levels = unique(tempx$Splice_Event))

    ##-- tile plot ---
    basesize <- 10
    ppx <- ggplot(data = tempx, aes(x=CI_Type, y=Splice_Event, fill=HR)) + 
    geom_tile()+
    scale_x_discrete()+
    scale_y_discrete()+
    xlab("Survival association type")+ylab("Splicing event")+
    theme_bw()+theme(axis.text.x = element_text(size = 0.6*basesize, angle = 60, vjust=0.5, hjust=0.5, colour = "black"),
        axis.text.y = element_text(size = 0.6*basesize, angle = 0, colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.title=element_text(size=basesize * 0.6, face="bold"),
        legend.title=element_text(size=basesize * 0.6), legend.text=element_text(size=basesize * 0.6))+coord_fixed()+
    guides(fill=guide_legend(title="Hazard ratio (log2)",ncol=1))
    ggsave(ppx,filename=paste0(save_dir, "/",all_cancer[k],".png"),width=4, height=3, dpi=400)

}



##--- some KM plots ------------------------------



