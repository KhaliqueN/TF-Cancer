##############################################################################################
# Purpose: Dependency map of the selected TFs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

save_dir <- '../results_new/Dependency'

if(!dir.exists(save_dir)){
	dir.create(save_dir, recursive=TRUE)
}

##-- TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
dep_map <- as.data.frame(data.table::fread('../data/CRISPR_(DepMap_Public_24Q4+Score,_Chronos)_subsetted.csv'))
all_lin2 <- unique(dep_map$lineage_3)
cols <- colnames(dep_map)
c2consi <- list("Bladder Urothelial",c("Breast Invasive","Invasive Breast"),"Colorectal Adenocarcinoma",
    "Head and Neck Squamous Cell Carcinoma",
    "Renal Clear Cell Carcinoma", "Hepatocellular Carcinoma", "Lung Adenocarcinoma", "Lung Squamous Cell Carcinoma", 
    "Prostate Adenocarcinoma",
    "Stomach Adenocarcinoma", "Thyroid Cancer", "Endometrial Carcinoma")
filt_cancer <- c('BLCA', 'BRCA','COAD','HNSC','KIRC', 'LIHC', 'LUAD','LUSC','PRAD','STAD','THCA','UCEC')
##----------------------------------------------------------------------

input_dir <- '../data/PSI_data'
fdr <- 0.05
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))

all_cancer <- substr(basename(all_files), 1,4)

##------------------------------------------------
tf_events <- list()
for(k in 1:length(filt_cancer)){

    temp <- data.table::fread(all_files[which(all_files %like% filt_cancer[k])], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy <- tempy[order(-abs(tempy$MEDIAN_DIFF)), ]
    tempz <- tempy[tempy$MEDIAN_DIFF > 0.2,] 
    tf_events[[k]] <- unique(tempz$symbol) 
}


##--- dependecy plots -----------------------------
color <- colorRampPalette(c("navy", "white", "red"))(13)
breaks <- seq(-1.5, 1.5, 0.25)

for(k in 1:length(filt_cancer)){

    # c("cell_line_display_name", which(cols %in% tf_events[[k]]))
    temp <- dep_map[which(dep_map$lineage_3 %like% c2consi[[k]]),]
    tempx <- temp[, which(cols %in% c("cell_line_display_name",tf_events[[k]]))]

    if(nrow(tempx) < 2){next}
    rownames(tempx) <- tempx$cell_line_display_name
    tempx <- tempx[,-1]
    col_names <- names(sort(apply(tempx, 1, function(x) length(which(x < 0)))))
    row_names <- names(sort(apply(tempx, 2, function(x) length(which(x < 0)))))

    tempx$colid <- rownames(tempx)
    library(dplyr)
    tempy <- tempx %>% arrange(colid, col_names)
    tempz <- tempy[,-length(tempy)]
    tempzz <- tempz[, match(rev(row_names), colnames(tempz))]

    pdx <- t(as.matrix(tempzz))
    p <- pheatmap(pdx, color=rev(color), breaks=breaks, fontsize=3, cluster_rows=FALSE, cluster_cols=FALSE)
    ggsave(p,filename=paste0(save_dir, "/",filt_cancer[k],".png"),width=5, height=3, dpi=400)

}

