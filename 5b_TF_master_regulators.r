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

##--- cancer non-specific master regulators
tmas <- as.data.frame(readxl::read_excel('../data/sciadv.abf6123_tables_s1_to_s14.xlsx', 7))
tmas4 <- tmas[-1,]
colnames(tmas4) <- c('Cancer','Symbol')
## -----------------------------------------------------------------------

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

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
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

    # ##----- expressoin of genes ---
    # temp_filed <- data.table::fread(diffExprFiles[k])
    # temp_filex <- temp_filed[temp_filed$Ensembl_gene_id %in% tfs$Ensembl_Gene_ID, ]
    # wh1 <- which(temp_filex$log2FoldChange < -2)
    # wh2 <- which(temp_filex$log2FoldChange > 2)
    # temp_filey <- temp_filex[union(wh2, wh1), ]
    # temp_filez <- temp_filey[temp_filey$padj < 0.05, ]

    # temp_file <- as.data.frame(data.table::fread(exprFiles[k]))
    # last_flag <- substr(colnames(temp_file), 14,15)
    # n_col1 <- which(last_flag == '11') ## control samples
    # n_col2 <- which(last_flag == '01') ## cancer samples
    # temp_file <- temp_file[, n_col2]
    # temp_file$mean <- rowMeans(temp_file)
    # temp_file$min <- apply(temp_file, 1, FUN = min)
    # temp_file$max <- apply(temp_file, 1, FUN = max)

    # temp_file$Ensembl_gene_id <- temp_filed$Ensembl_gene_id
    # temp_filexx <- temp_file[temp_file$Ensembl_gene_id %in% tfs$Ensembl_Gene_ID, ]
    # temp_filexx <- temp_filexx[,seq(length(temp_filexx)-3,length(temp_filexx))]
    # temp_filexx <- temp_filexx[order(-temp_filexx$mean), ]

    # temp_fileyy <- merge(temp_filexx, emap, by='Ensembl_gene_id')
    # temp_fileyy <- temp_fileyy[order(-temp_fileyy$mean), ]
    # temp_fileyy <- temp_fileyy[which(temp_fileyy$HGNC_symbol != ''), ]
    # tempz <- temp_fileyy[seq(1,ceiling(0.05*length(temp_fileyy[[1]]))),]
    # ##-----------------------------

    # top5xpr[[k]] <- intersect(tempz$HGNC_symbol, temptfs)
    top5xpro[[k]] <- intersect(tmas3[tmas3$Cancer == all_cancer[k], ]$Symbol, temptfs)
    # bt95xpr[[k]] <- intersect(tmas4[tmas4$Cancer == all_cancer[k], ]$Symbol, temptfs)

    ovrlap <- intersect(tempmas$Symbol[1:79], temptfs)
    tempz <- tempy[tempy$symbol %in% ovrlap, ]
    tempz <- tempz[order(-abs(tempz$MEDIAN_DIFF)), ]
    tsym <- c(tsym, tempz$symbol[1])
    tsty <- c(tsty, tempz$splice_type[1])
    tpos <- c(tpos, tempz$POS[1])
    tneg <- c(tneg, tempz$NEG[1])
    tdiff <- c(tdiff, tempz$MEDIAN_DIFF[1])
    tsid <- c(tsid, tempz$as_id[1])
}

##--- plot the number of splicing events affecting TFs -------------------
pdata <- data.frame(cancer=all_cancer, count1=lengths(top5pr), count2=lengths(top5xpro))
pdata <- pdata[-4,]
pdata$count3 <- pdata$count1-pdata$count2
pdata <- pdata[,-2]
colnames(pdata) <- c('Cancer','Top 5%','Others')
pdata <- reshape2::melt(pdata)
pdata$count <- rep(lengths(top5pr)[-4],2)
p <- ggplot(pdata, aes(Cancer, value, fill=variable)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of TFs", limits=c(0,max(pdata$count))) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
geom_text(aes(label=value), position=position_stack(vjust=0.5), size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Expression rank",ncol=1))
ggsave(p,filename=paste0(save_dir,"/Master_regulators.png"),width=5, height=3, dpi=400)



###-------- table for top events ---------
toptable <- data.frame(CANCER=all_cancer, GENE=tsym, SPLICE_ID=tsid, SPLICE_TYPE=tsty, MEDIAN_DIFF=tdiff, POS=tpos, NEG=tneg, PAIRED=paired_sam[[2]])

##----------------------------------------



##-- Look into the splicing behaviour in the master TFs -------------------
cancer <- c()
mcr <- c()
mnl <- c()
diff <- c()
asid <- c()
gene <- c()
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Master_regulators.xlsx'))

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    tgenes <- top5pr[[k]]

    mcrx <- c()
    mnlx <- c()
    diffx <- c()
    asidx <- c()
    genex <- c()

    if(length(tgenes) != 0){
        for(j in 1:length(tgenes)){
            tempy <- tempx[tempx$symbol == tgenes[j],]
            mcrx <- c(mcrx, tempy$MEDIAN_CANCER)
            mnlx <- c(mnlx, tempy$MEDIAN_NORMAL)
            diffx <- c(diffx, tempy$MEDIAN_DIFF)
            asidx <- c(asidx, tempy$as_id)
            genex <- c(genex, rep(tgenes[j], length(tempy[[1]])))

            mcr <- c(mcr, tempy$MEDIAN_CANCER)
            mnl <- c(mnl, tempy$MEDIAN_NORMAL)
            diff <- c(diff, tempy$MEDIAN_DIFF)
            cancer <- c(cancer, rep(all_cancer[k], length(tempy[[1]])))
            asid <- c(asid, tempy$as_id)
            gene <- c(gene, rep(tgenes[j], length(tempy[[1]])))
        }

        ##-- save excel sheet ----
        pdatax <- data.frame(gene=genex, asid=asidx, median_cancer=mcrx, median_normal=mnlx, median_diff=diffx)
        pdatax <- pdatax[order(-abs(pdatax$median_diff)),]
        openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
        openxlsx::writeData(wb1, sheet = all_cancer[k], pdatax)
        openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Master_regulators.xlsx'), overwrite = T)
    }
}

pdata <- data.frame(cancer=cancer, gene=gene, asid=asid, median_cancer=mcr, median_normal=mnl, median_diff=diff)
p <- ggplot(pdata, aes(cancer, median_diff)) + 
geom_jitter(aes(color=cancer),size=0.4)+
geom_violin()+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="Median of PSI differences \namong paired samples", limits=c(-1,1)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# geom_text(aes(label=value), position=position_stack(vjust=0.5), size=3)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a',
    '#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color='none')
ggsave(p,filename=paste0(save_dir,"/Master_regulators_splicingEvents.png"),width=7, height=4, dpi=400)






