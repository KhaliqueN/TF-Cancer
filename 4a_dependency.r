##############################################################################################
# Purpose: Dependency map of the selected TFs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

save_dir <- '../results_new/Dependency'
if(dir.exists(save_dir)){
    unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

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

all_files <- all_files[-4]
all_files_raw <- all_files_raw[-4]

all_cancer <- substr(basename(all_files), 1,4)

paired_sam <- data.table::fread('../data/cancer_paired_samples.txt')
paired_sam <- paired_sam[-4]

##------------------------------------------------
tf_events <- list()
tf_events_all <- list()
tf_max_medians <- list()
tf_max_events <- list()
for(k in 1:length(filt_cancer)){

    temp <- data.table::fread(all_files[which(all_files %like% filt_cancer[k])], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy <- tempy[order(-abs(tempy$MEDIAN_DIFF)), ]
    tempz <- tempy[tempy$MEDIAN_DIFF > 0.2,] 
    tf_events[[k]] <- unique(tempz$symbol) 
    ugs <- unique(tempy$symbol) 
    tf_events_all[[k]] <- ugs

    gvls <- c()
    evls <- c()
    for(j in 1:length(ugs)){
        tempyx <- tempy[tempy$symbol == ugs[j], ]
        if(nrow(tempyx) > 1){
            wh <- which(abs(tempyx$MEDIAN_DIFF) == max(abs(tempyx$MEDIAN_DIFF)))
            mnv <- tempyx$MEDIAN_DIFF[wh[1]] ## taking one of them
            mne <- paste0(tempyx$as_id[wh[1]],'_',tempyx$splice_type[wh[1]])
        }else{
            mnv <- tempyx$MEDIAN_DIFF
            mne <- paste0(tempyx$as_id,'_',tempyx$splice_type)
        }
        gvls <- c(gvls, mnv)
        evls <- c(evls, mne)
    }

    tf_max_medians[[k]] <- gvls
    tf_max_events[[k]] <- evls

}


##--- dependecy plots -----------------------------
# color <- colorRampPalette(c("navy", "white", "red"))(13)
color <- c(rev(c('#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d')),'#ffffff','#c6dbef','#9ecae1')
# breaks <- c(seq(-2.75, 0, 0.5), 0, 0.25, 0.5, 0.75)
breaks <- seq(-3, 1, 0.5)
amin <- c()
amax <- c()

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

    ##--- add splicing details -------------
    tempcol <- colnames(tempzz)
    ids <- c()
    for(j in 1:length(tempcol)){
        wh <- which(tf_events_all[[k]] == tempcol[j])
        ids <- c(ids, paste0(tf_events_all[[k]][wh], '_', tf_max_events[[k]][wh],' (',tf_max_medians[[k]][wh],')'))
    }

    colnames(tempzz) <- ids
    ##--------------------------------------

    amin <- c(amin, min(tempzz, na.rm=TRUE))
    amax <- c(amax, max(tempzz, na.rm=TRUE))

    pdx <- t(as.matrix(tempzz))
    p <- pheatmap(pdx, color=color, breaks=breaks, fontsize=3, cluster_rows=FALSE, cluster_cols=FALSE,cellheight=5, cellwidth = 5)
    ggsave(p,filename=paste0(save_dir, "/",filt_cancer[k],".png"),width=7, height=5, dpi=600)

}




all_sig_events <- unique(unlist(tf_events_all))
cvals <- rep(100,length(all_sig_events))
svals <- rep(100,length(all_sig_events))

for(k in 1:length(filt_cancer)){
    temp <- dep_map[which(dep_map$lineage_3 %like% c2consi[[k]]),]
    tempx <- temp[, which(cols %in% c("cell_line_display_name",tf_events_all[[k]]))]
    
    if(nrow(tempx) < 2){next}
    rownames(tempx) <- tempx$cell_line_display_name
    tempx <- tempx[,-1]
    acols <- colnames(tempx)

    for(j in 1:length(acols)){
        tx <- min(tempx[[j]], na.rm=TRUE)
        whp <- which(all_sig_events == acols[j])
        if(cvals[whp] > tx){
            cvals[whp] <- tx
        }
        wh <- which(tf_events_all[[k]] == acols[j])
        svals[whp] <- tf_max_medians[[k]][wh]
    }
}


pdata <- data.frame(GENE=all_sig_events, MINCHRONOS=cvals, MAXMEDIANDIFF=svals)
pdata <- pdata[pdata$MINCHRONOS != 100, ] ## 33 had no values in the chronos evaluations --> removing those
pdata <- pdata[pdata$MAXMEDIANDIFF != 100, ] 

basesize <- 10
ppx <- ggplot(data = pdata, aes(x=MINCHRONOS, y=MAXMEDIANDIFF)) + 
geom_point(alpha=0.6)+
scale_x_continuous()+
scale_y_continuous()+
geom_hline(yintercept=0.2, color='red', linetype='dashed')+
geom_hline(yintercept=-0.2, color='red', linetype='dashed')+
geom_vline(xintercept=0, color='blue', linetype='dashed')+
xlab("Chronos score")+ylab("Difference between \nmedians of paired samples")+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill='none')
ggsave(ppx,filename=paste0(save_dir, "/Scatter_dependency.png"),width=4, height=3, dpi=400)


## cumulative counts --
# pdata <- data.frame(GENE=all_sig_events, MINCHRONOS=cvals)
# pdata <- pdata[pdata$MINCHRONOS != 100, ] ## 33 had no values in the chronos evaluations --> removing those
# pdata$wbreak <- cut(pdata$MINCHRONOS, breaks = seq(-3.25,0.25,0.25), include.lowest=TRUE)

# pdatat <- plyr::count(pdata$wbreak)
# pdatat$cumsum <- cumsum(pdatat[[2]])

# basesize <- 10
# ppx <- ggplot(data = pdatat, aes(x=x, y=cumsum)) + 
# geom_bar(stat="identity")+
# scale_y_continuous(limits=c(0,max(pdatat$cumsum)+100), breaks = seq(0, max(pdatat$cumsum)+100, by = 100))+
# xlab("Chronos score range")+ylab("Cumulative # of TFs with \nat least one perturbed event")+
# geom_text(aes(label=cumsum), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
#     axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
#     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"))+
# guides(fill='none')
# ggsave(ppx,filename=paste0(save_dir, "/Cumulative_dependency.png"),width=4, height=3, dpi=400)

