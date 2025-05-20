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

if(!dir.exists(save_dir)){
	dir.create(save_dir, recursive=TRUE)
}

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

    # tempy$POSP <- tempy$POS/paired_sam[[2]][k]
    # tempy$NEGP <- tempy$NEG/paired_sam[[2]][k]
    # wh1 <- which(tempy$POSP > num_sams)
    # wh2 <- which(tempy$NEGP > num_sams)
    # wh <- union(wh1, wh2)
    # tempy <- tempy[wh,] 

    ## save excel sheet ----
    tempz <- data.frame(tempy)
    tdatat <- tempz[, c(1,2,3, seq(length(tempz)-11, length(tempz)))]
    tdatat <- tdatat[order(tdatat$FDR), ]

    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tdatat)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Perturbed_TF_splicing_events.xlsx'), overwrite = T)

    # ##--- bin the gained fraction ----
    # posp_break <- cut(tempy$POSP, breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1), include.lowest=TRUE)
    # negp_break <- cut(tempy$NEGP, breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1), include.lowest=TRUE)
    
    events_tf[[k]] <- tempy$as_id
    # txx <- plyr::count(posp_break)
    # txx$pos_per <- txx$freq/sum(txx$freq)
    # txx$prct <- signif((txx$freq/sum(txx$freq))*100, 3)

    # tyy <- plyr::count(negp_break)
    # tyy$neg_per <- tyy$freq/sum(tyy$freq)
    # tyy$prct <- signif((tyy$freq/sum(tyy$freq))*100, 3)

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

    # ##--- bin the gained fraction ----
    # posp_break <- cut(tempy$POSP, breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), include.lowest=TRUE)
    # negp_break <- cut(tempy$NEGP, breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), include.lowest=TRUE)
    
    # txx <- plyr::count(posp_break)
    # txx$pos_per <- txx$freq/sum(txx$freq)
    # txx$prct <- signif((txx$freq/sum(txx$freq))*100, 3)

    # tyy <- plyr::count(negp_break)
    # tyy$neg_per <- tyy$freq/sum(tyy$freq)
    # tyy$prct <- signif((tyy$freq/sum(tyy$freq))*100, 3)

}

# pdata <- data.frame(CANCER=tcancer, High=tpos, Low=tneg)
# pdatat <- reshape2::melt(pdata)
# p <- ggplot(pdatat, aes(CANCER, value, fill=variable)) + 
# geom_boxplot()+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Ccancer type") + 
# scale_y_continuous(name="Fraction of paired samples", limits=c(0,1)) +
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99',
#     '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# guides(fill=guide_legend(title="PSI in cancer \nvs. \nnornal sample",ncol=1))
# ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_dist.png"),width=5, height=3, dpi=400)


pdata <- data.frame(cancer=tcancer, median_diff=tdiff)
# pdata$pbreak <- cut(pdata$median_diff, breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), include.lowest=TRUE)
# tempx <- pdata[abs(pdata$median_diff) > 0.25, ]
# tempx <- plyr::count(tempx$cancer)
p <- ggplot(pdata, aes(cancer, median_diff)) + 
geom_jitter(aes(color=cancer),size=0.1)+
geom_violin(outlier.shape = NA)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="Median of PSI differences\n across paired samples", limits=c(-0.75,0.75)) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99',
    '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color='none')
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_dist.png"),width=5, height=3, dpi=400)


# pdata <- data.frame(cancer=tcancer, median_diff=tdiff)
# pdata$pbreak <- cut(pdata$median_diff, breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), include.lowest=TRUE)
# pdataxx <- plyr::count(pdata[,c(1,3)], c("cancer","pbreak"))
# ppx <- ggplot(data = pdataxx, aes(x=cancer, y=freq, fill=pbreak, group=cancer)) + 
# geom_bar(stat="identity")+
# scale_y_continuous(limits=c(0,max(lengths(events_tf))+50), breaks = seq(0, max(lengths(events_tf)), by = 100))+
# xlab("Cancer type")+ylab("% of significant AS events")+
# scale_fill_manual(values=c('#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4'))+
# theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
#     axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
#     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"))+
# guides(fill=guide_legend(title="Median PSI \ndifference"))
# ggsave(ppx,filename=paste0(save_dir, "/Sig_events_types_TFs_dist.png"),width=7, height=3, dpi=400)

###---------------------------------------------------------------------------------------------------------------------





##-------- Pairwise overlap of slicing events affecting TFs across cancer types --------------------------------------------
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
ggsave(p,filename=paste0(save_dir,"/Events_perturbing_overlap.png"),width=3.5, height=3, dpi=400)


##----- PSI profile of the pancancer (at least 12 cancer types) events -------------
events_to_consider <- num_events_combs[[12]]
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
pdx <- as.matrix(pdata_mat)

p <- pheatmap(pdx)
ggsave(p,filename=paste0(save_dir, "/Pancancer_events.png"),width=7, height=5, dpi=400)








# ##---- KEGG pathways analysis of the affected genes -------------------------------------------------
# ##---- Using TFLink ----
# tf_ensemb_map <- as.data.frame(data.table::fread('../data/TF_ensembl_uniprot.txt', sep='\t'))
# diff_expr_files <- gtools::mixedsort(list.files('../data/Diff_expr', full.names=TRUE))
# ensembl_gene_map <- data.table::fread('../data/ensembl_name_map.txt')
# ##---- grn without direction downloaded from here: https://tflink.net/download/
# grn <- as.data.frame(data.table::fread('../data/TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv'))
# grn_filt <- grn[grn$`Name.TF` %in% tfs$Gene_Symbol,]
# ##---- 14,439 edges ----

# all_kegg <- data.table::fread('../data/all_human_pathways.txt', sep='\t')
# colnames(all_kegg) <- c('Ensembl_gene_id','Uniprotswissprot','Entrez_id','Go_id','Pathways_name')
# all_path <- plyr::count(all_kegg$Go_id) ## 347 total KEGG pathways
# wh_path <- unique(all_path[which(all_path$freq > 2), ]$x)
# allkegg <- all_kegg[all_kegg$Go_id %in% wh_path, ] ### all pathways remain
# bgSize <- length(unique(allkegg$Uniprotswissprot))
# pval <- c()
# tcancer <- c()
# Kterm <- c()
# Ktermn <- c()
# num_gene <- c()

# events_to_consider <- num_events_combs[[10]]
# temp <- data.table::fread(all_files[1], sep='\t')
# tempg <- unique(temp[temp$as_id %in% events_to_consider, ]$symbol)

# temp_gene <- unique(tf_ensemb_map[tf_ensemb_map$HGNC_symbol %in% tempg, ]$Uniprotswissprot)
# temp_gene <- unique(grn_filt[grn_filt$`UniprotID.TF` %in% temp_gene, ]$`UniprotID.Target`)

# num_gene <- length(temp_gene)

# temp_kegg <- allkegg[allkegg$Uniprotswissprot %in% temp_gene, ]
# sampleSize <- length(unique(temp_kegg$Uniprotswissprot))
# tempid <- unique(temp_kegg$Go_id)

# for(j in 1:length(tempid)){
#     setA <- unique(allkegg[allkegg$Go_id == tempid[j], ]$Uniprotswissprot)
#     setB <- unique(temp_kegg[temp_kegg$Go_id == tempid[j], ]$Uniprotswissprot)
#     ## hypergeometric test
#     hyp <- phyper(length(setB)-1,length(setA),bgSize-length(setA), sampleSize,lower.tail = FALSE)
#     pval <- c(pval, hyp)
#     Kterm <- c(Kterm, tempid[j])
#     Ktermn <- c(Ktermn, unique(temp_kegg[temp_kegg$Go_id == tempid[j], ]$Pathways_name))
# }

# ktd <- data.frame(Kterm, Ktermn, pval)
# ktd$FDR <- p.adjust(ktd$pval, 'fdr')
# ktdx <- ktd[ktd$FDR < 0.05, ]

# #-- save excel file -----
# wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/KEGG_enrichment_10cancer.xlsx'))
# temp1a <- ktdx[order(ktdx$FDR), ]
# colnames(temp1a) <- c('KEGG ID', 'KEGG pathway name', 'pvalue', 'FDR')
# openxlsx::addWorksheet(wb1, sheetName = 'Atleast 10 cancer')
# openxlsx::writeData(wb1, sheet = 'Atleast 10 cancer', temp1a)
# openxlsx::saveWorkbook(wb1, paste0(save_dir,'/KEGG_enrichment_10cancer.xlsx'), overwrite = T)
# ##-----------------------







##******* Redo using COSMIC data *******************************************
# ##--- Overlap with oncogenes and cancer suppressor genes ----------------
# ##-- speed up --read all files in memory!!
# all_files_store <- list()
# for(k in 1:length(all_files)){
#     temp <- as.data.frame(data.table::fread(all_files[k], sep='\t'))
#     ll <- length(temp) - 3
#     all_files_store[[k]] <- read.table(all_files[k], colClasses = c('character','character','character',rep('NULL', ll)), header = TRUE, fill=TRUE)
# }

# gene_ovr <- list()
# for(k in 1:length(all_cancer)){
#     events_to_consider1 <- num_events_combs_pos[[k]]
#     events_to_consider2 <- num_events_combs_neg[[k]]
#     g_ovr <- c()
#     for(j in 1:length(all_cancer)){
#         # temp <- data.table::fread(all_files[k], sep='\t')
#         temp <- all_files_store[[j]]
#         temp1 <- temp[temp$as_id %in% events_to_consider1, ]
#         temp2 <- temp[temp$as_id %in% events_to_consider2, ]
#         g_ovr <- union(union(g_ovr, temp1$symbol), temp2$symbol)
#     }  
#     gene_ovr[[k]] <- g_ovr
#     cat('Combination',k,'of',length(all_cancer),'done\n')
# }


# ##---download tumor suppressor genes from: https://bioinfo.uth.edu/TSGene/download.cgi
# #--https://bioinfo.uth.edu/TSGene/Human_TSGs.txt ---
# xx <- "Human_TSGs.txt"
# system(paste0("wget -O ../data/",xx," https://bioinfo.uth.edu/TSGene/",xx))
# tsgs <- data.table::fread(paste0("../data/",xx))

# cancergenes <- data.table::fread('../data/cancerGeneList.tsv')
# ## keep genes with at least two evidence
# wh <- which(cancergenes$`# of occurrence within resources (Column J-P)` > 2)
# cancergenes <- cancergenes[wh, ]
# oncogenes <- cancergenes[cancergenes$`Is Oncogene` == 'Yes', ]
# supgenes <- cancergenes[cancergenes$`Is Tumor Suppressor Gene` == 'Yes', ]
# tf_onco <- vector(mode = "list", length = length(all_cancer))
# tf_sup <- vector(mode = "list", length = length(all_cancer))

# for(k in 1:length(all_cancer)){ ## overlap with the TFs in multiple cancer types
#     # tf_onco[[k]] <- intersect(oncogenes$`Hugo Symbol`, gene_ovr[[k]])
#     tf_sup[[k]] <- intersect(tsgs$GeneSymbol, gene_ovr[[k]])
# }


##----- unique splicing events -----------------------------------------------------
num_events_unq <- list()
num_events_unq_counts <- c()
unique_events_tf <- list()

for(k in 1:length(all_cancer)){
    temp_unq <- TF_splicing_events[[k]]
    temp_un <- c()
    for(j in 1:length(all_cancer)){
        if(j != k){
            temp_un <- union(temp_un, TF_splicing_events[[j]])
        }
    }
    temp_unq <- setdiff(temp_unq, temp_un)
    num_events_unq[[k]] <- temp_unq
    num_events_unq_counts <- c(num_events_unq_counts, length(temp_unq))
}

pdata <- data.frame(cancer=all_cancer, count=num_events_unq_counts)

p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(pdata$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of significant splicing events", limits=c(0,maxv+20)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"),
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),  
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Unique_sig_events_TFs.png"),width=3.5, height=3, dpi=400)


##----- PSI profile of the top ranking cancer-specific unique events -------------
# events_to_consider <- num_events_combs[[12]]
events_to_consider <- c()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    temp1 <- temp[temp$as_id %in% num_events_unq[[k]], ]
    wh <- which(temp1$MEDIAN_DIFF == max(temp1$MEDIAN_DIFF))
    events_to_consider <- c(events_to_consider, temp1$as_id[wh[1]])
}

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
pdx <- as.matrix(pdata_mat)

p <- pheatmap(pdx)
ggsave(p,filename=paste0(save_dir, "/Cancerspecific_events.png"),width=7, height=5, dpi=400)







# ##--- Overlap with the survival associated events ---------------
# ## Study: OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers, NAR, 2022
# ## Downloaded http://47.98.127.64:8080/beta/download?fileName=spliceseq_clinical_as_survival.csv.gz
# splicing <- data.table::fread('../data/spliceseq_clinical_as_survival.csv')
# splicing <- as.data.frame(splicing)
# events_tf_surv <- list()
# events_tf_surv_c <- c()
# events_tf_surv_p <- c()
# wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Perturbed_TF_splicing_events_survival.xlsx'))

# for(k in 1:length(all_cancer)){

#     events_tfx <- events_tf[[k]]
#     temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
#     temp1$SE <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
#     whi <- intersect(temp1$SE, events_tfx)
#     tempy <- temp1[temp1$SE %in% whi, ]

#     ##--
#     tempsi <- data.table::fread(all_files[k])
#     tempsi <- tempsi[tempsi$as_id %in% whi, ]

#     tmdif <- c()
#     tfdr <- c()

#     for(j in 1:length(tempy[[1]])){
#     	tmdif <- c(tmdif, tempsi$MEDIAN_DIFF[which(tempsi$as_id == tempy$SE[j])])
#     	tfdr <- c(tfdr, tempsi$FDR[which(tempsi$as_id == tempy$SE[j])])
#     }

#     tempy$MEDIAN_DIFF <- tmdif
#     tempy$FDR <- tfdr
#     tempy <- tempy[order(-abs(tempy$MEDIAN_DIFF)),]

#     events_tf_surv[[k]] <- intersect(tempy$SE, events_tfx)
#     events_tf_surv_c <- c(events_tf_surv_c, length(events_tf_surv[[k]]))
#     events_tf_surv_p <- c(events_tf_surv_p, (length(events_tf_surv[[k]])/length(events_tfx)*100))

#     openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
#     openxlsx::writeData(wb1, sheet = all_cancer[k], tempy)
#     openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Perturbed_TF_splicing_events_survival.xlsx'), overwrite = T)

# }

# pdata <- data.frame(cancer=all_cancer, count=events_tf_surv_c, pcount=events_tf_surv_p)
# p <- ggplot(pdata, aes(cancer, pcount)) + 
# geom_bar(stat="identity",position=position_dodge())+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="% of perturbed events affecting TFs\n that are also survival associated", limits=c(0,75)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
# ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_surv.png"),width=3.5, height=3, dpi=400)




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


# ##-- overlap with tissue-specific TFs ---
# tissue_specific <- data.table::fread('../data/main.txt')
# tis <- tissue_specific[tissue_specific$`Gene Type` == 'TF', ]
# tis <- tis[tis$`Cell Type` == 'Cancer cell', ]
# all_tissue <- unique(tis$`Tissue Type`)

# sel_tissue <- list('Bladder', 'Breast', 'Colon', 'Esophagus', c('Nasopharyngeal carcinoma','Salivary gland','Oral cavity'),
#     'Kidney', 'Kidney', 'Kidney', 'Liver', 'Lung', 'Lung', 'Prostate', 'Stomach', 'Thyroid', 'Uterus')
# seltis <- list()
# fract <- c()
# for(k in 1:length(all_cancer)){
#     temptis <- unique(tis[tis$`Tissue Type` %in% sel_tissue[[k]], ]$`Gene Name`)
#     seltis[[k]] <- intersect(tf_events[[k]], temptis)
# }


# sel_tissue <- list('Bladder', 'Breast', 'Colon', 'Kidney', 'Liver', 'Lung', 'Prostate', 'Stomach', 'Thyroid', 'Uterus')
# all_tfs <- c()
# cangrp <- list('BLCA','BRCA','COAD',c('KICH','KIRC','KIRP'),'LIHC',c('LUAD','LUSC'), 'PRAD', 'STAD','THCA','UCEC')
# for(k in 1:length(cangrp)){
#     tx <- c()
#     for(j in 1:length(cangrp[[k]])){
#         tx <- union(tx, tf_events[[which(all_cancer == cangrp[[k]][j])]])
#     }
#     all_tfs[[k]] <- tx
# }

# seltis <- list()
# fract <- c()
# for(k in 1:length(sel_tissue)){
#     temptis <- unique(tis[tis$`Tissue Type` %in% sel_tissue[[k]], ]$`Gene Name`)
#     seltis[[k]] <- intersect(all_tfs[[k]], temptis)
#     fract <- c(fract, length(intersect(all_tfs[[k]], temptis))/length(temptis))
# }


