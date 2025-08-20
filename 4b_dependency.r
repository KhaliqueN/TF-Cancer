##############################################################################################
# Purpose: Dependency map of the selected TFs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)
library(cowplot)

save_dir <- '../results_new/Dependency'
if(dir.exists(save_dir)){
    unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

##-- TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
## Downloaded from : https://depmap.org/portal/data_page/?tab=customDownloads
## CRISPR -- exclude NA -- all celllines -- all genes --add celline metadata -- merge into single file
dep_map <- as.data.frame(data.table::fread('../data/CRISPR_(DepMap_Public_25Q2+Score,_Chronos)_subsetted.csv'))
all_lin2 <- unique(dep_map$lineage_3)
cols <- colnames(dep_map)
c2consi <- list("Bladder Urothelial","Breast Invasive Carcinoma","Colorectal Adenocarcinoma",
    "Head and Neck Squamous Cell Carcinoma",
    "Renal Clear Cell Carcinoma", "Hepatocellular Carcinoma", "Lung Adenocarcinoma", "Lung Squamous Cell Carcinoma", 
    "Prostate Adenocarcinoma","Stomach Adenocarcinoma", "Thyroid Cancer", "Endometrial Carcinoma")
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
##------------------------------------------------
# tf_events <- list()
tf_events_all <- list()
tf_max_medians <- list()
tf_max_events <- list()

for(k in 1:length(filt_cancer)){

    temp <- data.table::fread(all_files[which(all_files %like% filt_cancer[k])], sep='\t')
    wha <- which(temp$FDR < fdr)
    whb <- which(abs(temp$MEAN_NORMAL-temp$MEAN_CANCER) > fdr)
    wh <- intersect(wha, whb)
    tempx <- temp[wh, ]

    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy <- tempy[order(-abs(tempy$MEAN_DIFF)), ]
    # tempz <- tempy[tempy$MEAN_DIFF > 0.2,] 
    # tf_events[[k]] <- unique(tempz$symbol) 
    ugs <- unique(tempy$symbol) 
    tf_events_all[[k]] <- ugs

    gvls <- c()
    evls <- c()
    for(j in 1:length(ugs)){
        tempyx <- tempy[tempy$symbol == ugs[j], ]
        if(nrow(tempyx) > 1){
            wh <- which(abs(tempyx$MEAN_DIFF) == max(abs(tempyx$MEAN_DIFF)))
            mnv <- tempyx$MEAN_DIFF[wh[1]] ## taking one of them
            mne <- paste0(tempyx$as_id[wh[1]],'_',tempyx$splice_type[wh[1]])
        }else{
            mnv <- tempyx$MEAN_DIFF
            mne <- paste0(tempyx$as_id,'_',tempyx$splice_type)
        }
        gvls <- c(gvls, mnv)
        evls <- c(evls, mne)
    }

    tf_max_medians[[k]] <- gvls
    tf_max_events[[k]] <- evls

}


all_sig_events <- unique(unlist(tf_events_all)) ## 688 TFs related to all PTSEs in the 12 cancer types
pdata <- data.frame(matrix(ncol=5, nrow=0))
for(k in 1:length(filt_cancer)){
    temp <- dep_map[which(dep_map$lineage_3 %like% c2consi[[k]]),]
    tempx <- temp[, which(cols %in% c("cell_line_display_name",tf_events_all[[k]]))]
    
    if(nrow(tempx) < 1){next}
    rownames(tempx) <- tempx$cell_line_display_name
    tempx <- tempx[,-1]
    acols <- colnames(tempx)

    cvals <- rep(100,length(all_sig_events))
    svals <- rep(100,length(all_sig_events))
    celline <- rep(100,length(all_sig_events))
    tcancer <- rep(100,length(all_sig_events))
    asid <- rep(100,length(all_sig_events))

    for(j in 1:length(acols)){
        nac <- length(which(is.na(tempx[[j]])))
        if(nac == length(tempx[[j]])){
            next
        }else if(length(tempx[[j]]) == 1){
            tx <- tempx[[j]]
        }else{
            tx <- min(tempx[[j]], na.rm=TRUE)
        }
        whm <- which(tempx[[j]] == tx)
        whp <- which(all_sig_events == acols[j])
        cvals[whp] <- tx
        wh <- which(tf_events_all[[k]] == acols[j])
        svals[whp] <- tf_max_medians[[k]][wh]
        celline[whp] <- paste(rownames(tempx)[whm],collapse=',')
        tcancer[whp] <- filt_cancer[k]
        asid[whp] <- strsplit(tf_max_events[[k]][wh],'[_]')[[1]][1]
    }

    pdatat <- data.frame(CANCER=tcancer, GENE=all_sig_events, MINCHRONOS=cvals, MAXMEANDIFF=svals, CELLLINE=celline, ASID=asid)
    pdatat <- pdatat[pdatat$MINCHRONOS != 100, ] 
    pdata <- rbind(pdata, pdatat)

}


# length(which(pdata$MINCHRONOS == 100))
# pdata <- pdata[pdata$MAXMEANDIFF != 100, ] 
## 32 TFs had no values in the chronos evaluations --> removing those --> 669 TFs remain

pdata$col <- ifelse(pdata$MAXMEANDIFF > 0, 'A', 'B')
basesize <- 10
ppx <- ggplot(data = pdata, aes(y=MINCHRONOS, x=MAXMEANDIFF, color=col)) + 
geom_point(size=0.75, alpha=0.75)+
scale_x_continuous(limits=c(-0.8,0.8))+
scale_y_continuous()+
# geom_hline(yintercept=0.2, color='red', linetype='dashed')+
# geom_hline(yintercept=-0.2, color='red', linetype='dashed')+
# geom_vline(xintercept=0, color='blue', linetype='dashed')+
scale_color_manual(values=c('#d95f02','#1b9e77'))+
ylab("Chronos score")+xlab("Mean \u0394PSI")+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(color='none')
ggsave(ppx,filename=paste0(save_dir, "/Scatter_dependency.png"),width=2, height=2.5, dpi=500)

# pdatax <- pdata[abs(pdata$MAXMEDIANDIFF) > 0.2,  ]
# pdatax <- pdatax[order(pdatax$MINCHRONOS), ]

##--- Fraction of patients affected ------
max_chronos <- -0.5
pdataxx <- data.frame(matrix(ncol=7, nrow=0))
for(k in 1:length(filt_cancer)){

    temp <- data.table::fread(all_files[which(all_files %like% filt_cancer[k])], sep='\t')
    wha <- which(temp$FDR < fdr)
    whb <- which(abs(temp$MEAN_NORMAL-temp$MEAN_CANCER) > fdr)
    wh <- intersect(wha, whb)
    tempx <- temp[wh, ]
    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]

    pdatay <- pdata[pdata$CANCER == filt_cancer[k], ]
    pdatay <- pdatay[pdatay$MINCHRONOS < max_chronos, ]

    if(nrow(pdatay) == 0){next}

    ## add fraction of patients affected --
    whp <- which(paired_sam[[1]] == filt_cancer[k])
    sams <- c()
    stype <- c()
    for(j in 1:length(pdatay[[1]])){
        sams <- c(sams, max(tempy[tempy$as_id == pdatay$ASID[j], ]$POS, tempy[tempy$as_id == pdatay$ASID[j], ]$NEG))
        stype <- c(stype, tempy[tempy$as_id == pdatay$ASID[j], ]$splice_type)
    }

    pdatay$PAT_FRAC <- (sams/paired_sam[[2]][whp])*100
    pdatay$STYPE <- stype
    pdataxx <- rbind(pdataxx, pdatay)
}

pdataxx$ID <- paste0(pdataxx$GENE,'_',pdataxx$ASID,'_',pdataxx$STYPE,'\n','(',pdataxx$CANCER,', ',signif(pdataxx$MAXMEANDIFF,3),')')

tolabel <- subset(pdataxx, MINCHRONOS < -2)$ID
whl <- which(pdataxx$ID %in% tolabel)
pdataxx$PL <- ""
pdataxx$PL[whl] <- tolabel

p <- ggplot(pdataxx, aes(PAT_FRAC, MINCHRONOS, color=CANCER, label=PL)) + 
geom_point()+
theme(legend.text=element_text(size=12))
basesize <- 10
p <- p + 
scale_x_continuous(name="% of patients", limits=c(0,110), breaks = seq(0, 100, by = 20)) + 
scale_y_continuous(name="Minimum Chronos score") +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#fb9a99',
    '#fdbf6f','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
geom_text_repel(family = "Poppins",
    max.overlaps=Inf,
                      size = 2,
                      color='black',
                      arrow = arrow(length = unit(0.010, "npc")),
                      min.segment.length = 0) +
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"), legend.position="bottom",axis.title=element_text(size=basesize),
    legend.text=element_text(size=basesize), legend.title=element_text(size=basesize),
    panel.border = element_blank())+
guides(color='none')
ggsave(p,filename=paste0(save_dir,"/Chronos_vs_patfrac.png"),width=4, height=3, dpi=500)


##--- dependency plots -----------------------------
# color <- colorRampPalette(c("navy", "white", "red"))(13)
color <- c(rev(c('#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d')),'#ffffff','#c6dbef','#9ecae1')
# breaks <- c(seq(-2.75, 0, 0.5), 0, 0.25, 0.5, 0.75)
breaks <- seq(-3, 1, 0.5)
amin <- c()
amax <- c()
ncells <- c()
ntfs <- c()
all_plots <- list()

for(k in 1:length(filt_cancer)){

    whhc <- c()
    for(j in 1:length(c2consi[[k]])){
        whhc <- which(dep_map$lineage_3 %like% c2consi[[k]][j])
    }

    temp <- dep_map[whhc,]
    ncells <- c(ncells, length(temp[[1]]))
    whhc <- which(cols %in% c("cell_line_display_name",tf_events_all[[k]]))
    ntfs <- c(ntfs, length(whhc)-1)

    if(length(whhc) < 2){next}

    tempx <- temp[, whhc]

    if(nrow(tempx) < 1){next}
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
    temp_val <- c()
    for(j in 1:length(tempcol)){
        wh <- which(tf_events_all[[k]] == tempcol[j])
        ids <- c(ids, paste0(tf_events_all[[k]][wh], '_', tf_max_events[[k]][wh],' (',signif(tf_max_medians[[k]][wh],3),')'))
        temp_val <- c(temp_val, signif(tf_max_medians[[k]][wh],3))
    }

    colnames(tempzz) <- ids

    ## top 1% --------------
    names(temp_val) <- ids
    temp_val <- rev(sort(abs(temp_val)))
    
    cwh <- round(length(temp_val)*0.1)
    wh <- which(colnames(tempzz) %in% names(temp_val[1:cwh]))
    tempzz <- tempzz[,wh]
    ##--------------------------------------
    ## remove ptses with positive chronos scores
    lcs <- sapply(tempzz, function(x) which(x < 0))
    lcsn <- names(which(lengths(lcs) == 0))
    whc <- which(colnames(tempzz) %in% lcsn)
    if(length(whc) > 0){tempzz <- tempzz[,-whc]}
    
    amin <- c(amin, min(tempzz, na.rm=TRUE))
    amax <- c(amax, max(tempzz, na.rm=TRUE))

    pdx <- t(as.matrix(tempzz))
    p1 <- pheatmap(pdx, color=color, breaks=breaks, fontsize=4, cluster_rows=FALSE, cluster_cols=FALSE,cellheight=5, cellwidth = 5)
    # all_plots[[k]] <- p1[[4]] 
    ggsave(p1[[4]],filename=paste0(save_dir, "/",filt_cancer[k],".png"),width=6, height=6, dpi=600)

}

# combined_plot <- cowplot::plot_grid(plotlist=all_plots,labels=filt_cancer) #nrow & ncol depend on how you want to 
# ggsave(combined_plot,filename=paste0(save_dir, "/Dependency.png"),width=7, height=8, dpi=600)
         
                          

## dependecy all ----
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/PTSE_dependency.xlsx'))

for(k in 1:length(filt_cancer)){

    temp <- pdata[pdata$CANCER == filt_cancer[k], ]
    if(nrow(temp) == 0){next}
    asid <- c()
    astype <- c()
    tempcol <- temp$GENE
    for(j in 1:length(tempcol)){
        wh <- which(tf_events_all[[k]] == tempcol[j])
        ids <- c(ids, paste0(tf_events_all[[k]][wh], '_', tf_max_events[[k]][wh],' (',signif(tf_max_medians[[k]][wh],3),')'))
        asid <- c(asid, strsplit(tf_max_events[[k]][wh],'[_]')[[1]][1])
        astype <- c(astype, strsplit(tf_max_events[[k]][wh],'[_]')[[1]][2])
    }

    temp$asid <- asid
    temp$astype <- astype
    temp <- temp[,-6]
    colnames(temp) <- c('CANCER', 'TF','MIN_CHRONOS','MAX_MEAN_DIFF','CELL_LINE','AS_ID','AS_TYPE')
    openxlsx::addWorksheet(wb1, sheetName = filt_cancer[k])
    openxlsx::writeData(wb1, sheet = filt_cancer[k], temp)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/PTSE_dependency.xlsx'), overwrite = T)
}

