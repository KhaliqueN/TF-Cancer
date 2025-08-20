##############################################################################################
# Purpose: Overlap of significant TF splicing events with cancer type-related master TFs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)

save_dir <- '../results_new/MRs'
if(dir.exists(save_dir)){
    unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

##-- TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
##----------------------------------------------------------------------

paired_sam <- data.table::fread('../data/cancer_paired_samples.txt')
paired_sam <- paired_sam[-4]

input_dir <- '../data/PSI_data'
fdr <- 0.05
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))
all_files <- all_files[-4]
all_files_raw <- all_files_raw[-4]
all_cancer <- substr(basename(all_files), 1,4)

## Master regulator data ---------------
## From this study: Predicting master transcription factors from pan-cancer expression data ----
## Higher CaCTS score is better ---
tmas <- as.data.frame(readxl::read_excel('../data/sciadv.abf6123_tables_s1_to_s14.xlsx', 3))
coln <- as.character(tmas[2,])
wh <- which(coln %in% all_cancer)
tmas1 <- tmas[,c(1,wh)]
coln <- as.character(tmas1[2,])
wh <- which(coln %in% all_cancer)
coln <- coln[wh]
colnames(tmas1) <- c('Symbol',coln)
tmas2 <- tmas1[-c(1,2,3),]

# cancer-specific master regulators
tmas <- as.data.frame(readxl::read_excel('../data/sciadv.abf6123_tables_s1_to_s14.xlsx', 6))
tmas <- tmas[-c(1,2),]
tmas3 <- tmas[,c(1,2,3)]
colnames(tmas3) <- c('Cancer','Symbol','Cacts')

# ##--- cancer non-specific master regulators
# tmas <- as.data.frame(readxl::read_excel('../data/sciadv.abf6123_tables_s1_to_s14.xlsx', 7))
# tmas4 <- tmas[-1,]
# colnames(tmas4) <- c('Cancer','Symbol')
# ## -----------------------------------------------------------------------

# ##-- Mean expression of TFs according to RNAseq data downloaded by me ---
# save_dirx <- '../data/Diff_expr'
# save_dirf <- '../../../public_data/TCGA_FPKM_counts'
# tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
# emap <- data.table::fread('../data/ensembl_name_map.txt')
# emap <- unique(emap[,c(1,8)])
# diffExprFiles <- list.files(save_dirx, full.names=TRUE) 
# exprFiles <- list.files(save_dirf, full.names=TRUE) 
# ##------------------------------------------------------------------------


##--- overlap with TFs showing significant splicing events ---------------
top5pr <- list()
# top5xpr <- list()
top5xpro <- list()

tsym <- c()
tsty <- c()
tpos <- c()
tneg <- c()
tdiff <- c()
tsid <- c()
tcancer <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wha <- which(temp$FDR < fdr)
    whb <- which(abs(temp$MEAN_NORMAL-temp$MEAN_CANCER) > fdr)
    wh <- intersect(wha, whb)
    tempx <- temp[wh, ]

    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    temptfs <- gtools::mixedsort(unique(tempy$symbol))
    
    wh <- which(colnames(tmas2) == all_cancer[k])
    if(length(wh) == 0){
        top5pr[[k]] <- ''
        # top5xpr[[k]] <- ''
        top5xpro[[k]] <- ''
        next
    }
    tempmas <- tmas2[,c(1,wh)]
    colnames(tempmas) <- c('Symbol','Cacts')
    tempmas <- tempmas[order(tempmas$Cacts, decreasing=TRUE),]
    # temp_median <- median(as.numeric(tempmas$Cacts))

    ##---------
    top5pr[[k]] <- intersect(tempmas$Symbol[1:79], temptfs) ## top 5% of the Cacts score

    # ##-----------------------------

    # top5xpr[[k]] <- intersect(tempz$HGNC_symbol, temptfs)
    top5xpro[[k]] <- intersect(tmas3[tmas3$Cancer == all_cancer[k], ]$Symbol, temptfs)
    # bt95xpr[[k]] <- intersect(tmas4[tmas4$Cancer == all_cancer[k], ]$Symbol, temptfs)

    ovrlap <- intersect(tmas3[tmas3$Cancer == all_cancer[k], ]$Symbol, temptfs)#intersect(tempmas$Symbol[1:79], temptfs)
    if(length(ovrlap) != 0){
        for(j in 1:length(ovrlap)){
            tempz <- tempy[tempy$symbol %in% ovrlap[j], ]
            tempz <- tempz[order(-abs(tempz$MEAN_DIFF)), ]
            tsym <- c(tsym, tempz$symbol[1])
            tsty <- c(tsty, tempz$splice_type[1])
            tpos <- c(tpos, tempz$POS[1])
            tneg <- c(tneg, tempz$NEG[1])
            tdiff <- c(tdiff, tempz$MEAN_DIFF[1])
            tsid <- c(tsid, tempz$as_id[1])
            tcancer <- c(tcancer, all_cancer[k])
        }
    }
}

##--- plot the number of splicing events affecting TFs -------------------
pdata <- data.frame(cancer=all_cancer, count1=lengths(top5pr), count2=lengths(top5xpro))
# pdata <- pdata[-4,]
pdata$count3 <- pdata$count1-pdata$count2
pdata <- pdata[,-2]
colnames(pdata) <- c('Cancer','Top 5%','Others')
pdata <- reshape2::melt(pdata)
pdata$count <- rep(lengths(top5pr),2)
pdata[pdata==0] <- NA
pdata$Gene <- c(unlist(lapply(top5xpro,function(x) paste(sort(x), collapse='\n'))), rep(NA, length(pdata[[1]])/2))

p <- ggplot(pdata, aes(Cancer, value, fill=variable)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 10
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of TFs", limits=c(0,max(pdata$count)+10)) +
# geom_text(aes(label=Gene), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
geom_text(aes(x=Cancer, y=count, label=Gene), position=position_stack(vjust=0), size=2)+
theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8), legend.position=c(0.7,0.8))+
guides(fill=guide_legend(title="Gene expression",ncol=1))
ggsave(p,filename=paste0(save_dir,"/Master_regulators.png"),width=5, height=3, dpi=500)



###-------- table for top events ---------
toptable <- data.frame(CANCER=tcancer, GENE=tsym, SPLICE_ID=tsid, SPLICE_TYPE=tsty, MEAN_DIFF=tdiff, POS=tpos, NEG=tneg)
toptable$event <- paste0(toptable$GENE,'_',toptable$SPLICE_ID,'_',toptable$SPLICE_TYPE)


pdata_mat <- as.data.frame(matrix(nrow=length(unique(toptable$event)), ncol=length(unique(toptable$CANCER)),NA))
rownames(pdata_mat) <- gtools::mixedsort(unique(toptable$event))
colnames(pdata_mat) <- gtools::mixedsort(unique(toptable$CANCER))

for(k in 1:length(toptable[[1]])){
    wh1 <- which(rownames(pdata_mat) == toptable$event[k])
    wh2 <- which(colnames(pdata_mat) == toptable$CANCER[k])
    pdata_mat[wh1, wh2] <- as.numeric(toptable$MEAN_DIFF[k])
}

# pdata_mat <- cbind(data.frame(Event=unique(pdata$event)), pdata_mat)
pdx <- as.matrix(pdata_mat)#t(as.matrix(pdata_mat))

p <- pheatmap(pdx,fontsize=3, cluster_rows=FALSE, cluster_cols=FALSE,cellheight=5, cellwidth = 5, border_color=NA)
ggsave(p,filename=paste0(save_dir, "/MR_events.png"),width=5, height=7, dpi=600)




###---- dependency plots -----------------------------

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

##------------------------------------------------
tf_events <- list()
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
    tempz <- tempy[tempy$MEAN_DIFF > 0.2,] 
    tf_events[[k]] <- unique(tempz$symbol) 
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
cvals <- rep(100,length(all_sig_events))
svals <- rep(100,length(all_sig_events))
celline <- rep(100,length(all_sig_events))
tcancer <- rep(100,length(all_sig_events))

for(k in 1:length(filt_cancer)){
    temp <- dep_map[which(dep_map$lineage_3 %like% c2consi[[k]]),]
    tempx <- temp[, which(cols %in% c("cell_line_display_name",tf_events_all[[k]]))]
    
    if(nrow(tempx) < 1){next}
    rownames(tempx) <- tempx$cell_line_display_name
    tempx <- tempx[,-1]
    acols <- colnames(tempx)

    for(j in 1:length(acols)){
        tx <- min(tempx[[j]], na.rm=TRUE)
        whm <- which(tempx[[j]] == tx)
        whp <- which(all_sig_events == acols[j])
        if(cvals[whp] > tx){
            cvals[whp] <- tx
            wh <- which(tf_events_all[[k]] == acols[j])
            svals[whp] <- tf_max_medians[[k]][wh]
            celline[whp] <- paste(rownames(tempx)[whm],collapse=',')
            tcancer[whp] <- filt_cancer[k]
        }
    }
}


pdata <- data.frame(CANCER=tcancer, GENE=all_sig_events, MINCHRONOS=cvals, MAXMEANDIFF=svals, CELLLINE=celline)
pdata <- pdata[pdata$MINCHRONOS != 100, ] 
pdata <- pdata[pdata$MAXMEANDIFF != 100, ] 
## 19 TFs had no values in the chronos evaluations --> removing those --> 669 TFs remain

pdata$col <- ifelse(pdata$MAXMEANDIFF > 0, 'A', 'B')
basesize <- 6
ppx <- ggplot(data = pdata, aes(x=MINCHRONOS, y=MAXMEANDIFF, color=col)) + 
geom_point(alpha=0.6, size=0.5)+
scale_x_continuous()+
scale_y_continuous()+
# geom_hline(yintercept=0.2, color='red', linetype='dashed')+
# geom_hline(yintercept=-0.2, color='red', linetype='dashed')+
# geom_vline(xintercept=0, color='blue', linetype='dashed')+
scale_color_manual(values=c('#d95f02','#1b9e77'))+
xlab("Chronos score")+ylab("Mean \u0394PSI")+coord_flip()+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(color='none')
ggsave(ppx,filename=paste0(save_dir, "/Scatter_dependency.png"),width=1.5, height=2.5, dpi=400)

# pdatax <- pdata[abs(pdata$MAXMEDIANDIFF) > 0.2,  ]
# pdatax <- pdatax[order(pdatax$MINCHRONOS), ]


##--- dependency plots -----------------------------
# color <- colorRampPalette(c("navy", "white", "red"))(13)
color <- c(rev(c('#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d')),'#ffffff','#c6dbef','#9ecae1')
# breaks <- c(seq(-2.75, 0, 0.5), 0, 0.25, 0.5, 0.75)
breaks <- seq(-3, 1, 0.5)
amin <- c()
amax <- c()

for(k in 1:length(filt_cancer)){

    whhc <- c()
    for(j in 1:length(c2consi[[k]])){
        whhc <- which(dep_map$lineage_3 %like% c2consi[[k]][j])
    }

    temp <- dep_map[whhc,]
    whhc <- which(cols %in% c("cell_line_display_name",tf_events_all[[k]]))

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

    ## top 5% --------------
    names(temp_val) <- ids
    temp_val <- rev(sort(abs(temp_val)))
    
    cwh <- round(length(temp_val)*0.1)
    wh <- which(colnames(tempzz) %in% names(temp_val[1:cwh]))
    tempzz <- tempzz[,wh]
    ##--------------------------------------

    amin <- c(amin, min(tempzz, na.rm=TRUE))
    amax <- c(amax, max(tempzz, na.rm=TRUE))

    pdx <- t(as.matrix(tempzz))
    p <- pheatmap(pdx, color=color, breaks=breaks, fontsize=4, cluster_rows=FALSE, cluster_cols=FALSE,cellheight=5, cellwidth = 5)
    ggsave(p,filename=paste0(save_dir, "/",filt_cancer[k],".png"),width=6, height=4, dpi=600)

}





# ##--- dependency plots MRs -----------------------------
# # color <- colorRampPalette(c("navy", "white", "red"))(13)
# color <- c(rev(c('#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d')),'#ffffff','#c6dbef','#9ecae1')
# # breaks <- c(seq(-2.75, 0, 0.5), 0, 0.25, 0.5, 0.75)
# breaks <- seq(-3, 1, 0.5)
# amin <- c()
# amax <- c()

# for(k in 1:length(filt_cancer)){

#     # c("cell_line_display_name", which(cols %in% tf_events[[k]]))
#     whhc <- c()
#     for(j in 1:length(c2consi[[k]])){
#         whhc <- which(dep_map$lineage_3 %like% c2consi[[k]][j])
#     }

#     temp <- dep_map[whhc,]
#     whhc <- which(cols %in% c("cell_line_display_name",top5xpro[[k]]))

#     if(length(whhc) < 2){next}

#     tempx <- temp[, whhc]
#     rownames(tempx) <- tempx$cell_line_display_name
#     tempx <- tempx[,-1]
#     col_names <- names(sort(apply(tempx, 1, function(x) length(which(x < 0)))))
#     row_names <- names(sort(apply(tempx, 2, function(x) length(which(x < 0)))))

#     tempx$colid <- rownames(tempx)
#     library(dplyr)
#     tempy <- tempx %>% arrange(colid, col_names)
#     tempz <- tempy[,-length(tempy)]
#     tempzz <- tempz[, match(rev(row_names), colnames(tempz))]

#     ##--- add splicing details -------------
#     tempcol <- colnames(tempzz)
#     ids <- c()
#     for(j in 1:length(tempcol)){
#         wh <- which(tf_events_all[[k]] == tempcol[j])
#         ids <- c(ids, paste0(tf_events_all[[k]][wh], '_', tf_max_events[[k]][wh],' (',signif(tf_max_medians[[k]][wh],3),')'))
#     }

#     colnames(tempzz) <- ids
#     ##--------------------------------------

#     amin <- c(amin, min(tempzz, na.rm=TRUE))
#     amax <- c(amax, max(tempzz, na.rm=TRUE))

#     pdx <- t(as.matrix(tempzz))
#     p <- pheatmap(pdx, color=color, breaks=breaks, fontsize=3, cluster_rows=FALSE, cluster_cols=FALSE,cellheight=5, cellwidth = 5)
#     ggsave(p,filename=paste0(save_dir, "/",filt_cancer[k],".png"),width=7, height=8, dpi=600)

# }







# ###-------- table for top events ---------
# tsym <- c()
# tsty <- c()
# tpos <- c()
# tneg <- c()
# tdiff <- c()
# tsid <- c()
# tcell <- c()
# tchr <- c()

# for(k in 1:length(filt)){

#     temp <- data.table::fread(all_files[k], sep='\t')
#     wh <- which(temp$FDR < fdr)
#     tempx <- temp[wh, ]
#     whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
#     tempy <- tempx[whx,]

#     ovrlap <- pdata[pdata$CANCER == filt_cancer[k], ]
#     ovrlap <- ovrlap[order(ovrlap$MINCHRONOS),]

#     tempz <- tempy[tempy$symbol %in% ovrlap, ]
#     tempz <- tempz[order(-abs(tempz$MEDIAN_DIFF)), ]
#     tsym <- c(tsym, tempz$symbol[1])
#     tsty <- c(tsty, tempz$splice_type[1])
#     tpos <- c(tpos, tempz$POS[1])
#     tneg <- c(tneg, tempz$NEG[1])
#     tdiff <- c(tdiff, tempz$MEDIAN_DIFF[1])
#     tsid <- c(tsid, tempz$as_id[1])
#     tchr <- c(tchr, )

# }

# ###---------------------------------------

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




