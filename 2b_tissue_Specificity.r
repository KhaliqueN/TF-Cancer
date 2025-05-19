##############################################################################################
# Purpose: Overlap with tissue-specific TFs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

save_dir <- '../results_new/TF_splicing'

if(!dir.exists(save_dir)){
	dir.create(save_dir, recursive=TRUE)
}

##-- TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
tissue_specific <- data.table::fread('../data/main.txt')
# paired_sam <- data.table::fread('../data/cancer_paired_samples.txt')
##----------------------------------------------------------------------

input_dir <- '../data/PSI_data'
fdr <- 0.05
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))

all_cancer <- substr(basename(all_files), 1,4)

##------------------------------------------------
##--- Number of TFs with perturbed splicing events ----------------------------------------------------------------------
tf_events <- list()
for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy <- tempy[order(-abs(tempy$MEDIAN_DIFF)), ]
    ##-- take top 5 % ---
    tempz <- tempy[tempy$MEDIAN_DIFF > 0.25,]
    tf_events[[k]] <- unique(tempz$symbol) 
}


##-- overlap with tissue-specific TFs ---
tissue_specific <- data.table::fread('../data/main.txt')
tis <- tissue_specific[tissue_specific$`Gene Type` == 'TF', ]
tis <- tis[tis$`Cell Type` == 'Cancer cell', ]
all_tissue <- unique(tis$`Tissue Type`)

sel_tissue <- list('Bladder', 'Breast', 'Colon', 'Esophagus', c('Nasopharyngeal carcinoma','Salivary gland','Oral cavity'),
    'Kidney', 'Kidney', 'Kidney', 'Liver', 'Lung', 'Lung', 'Prostate', 'Stomach', 'Thyroid', 'Uterus')
seltis <- list()
fract <- c()
for(k in 1:length(all_cancer)){
    temptis <- unique(tis[tis$`Tissue Type` %in% sel_tissue[[k]], ]$`Gene Name`)
    seltis[[k]] <- intersect(tf_events[[k]], temptis)
}

##-- plot --
temp_names <- c()
for(k in 1:length(seltis)){
    if(length(seltis[[k]]) > 1){
        temp_names <- c(temp_names, paste(seltis[[k]], collapse='\n'))
    }else{
        temp_names <- c(temp_names, '')
    }
}

pdata <- data.frame(cancer=all_cancer, count=temp_names)
pdata$fr <- rep(0,length(pdata[[1]]))
p <- ggplot(pdata, aes(cancer, fr)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 10
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="Tissue-specific TFs affected \nby perturbed splicing events", limits=c(0,70)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0.5, vjust=0, angle=0, size=2.1)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = 0, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Tissue_specific_TFs.png"),width=7, height=3.5, dpi=400)




##--- characteristics of the splicing events of these tissue-specific TFs ----------------------------
temp_events <- list()
for(k in 1:length(all_cancer)){

    if(length(seltis[[k]]) != 0){
        temp <- data.table::fread(all_files[k], sep='\t')
        wh <- which(temp$FDR < fdr)
        tempx <- temp[wh, ]
        whx <- which(tempx$symbol %in% seltis[[k]]) ## number of AS events concerning TFs
        tempy <- tempx[whx,]
        temp_events[[k]] <- tempy
    }else{
        temp_events[[k]] <- ''
    }

}










