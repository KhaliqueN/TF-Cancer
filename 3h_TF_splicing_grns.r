##############################################################################################
# Purpose: Evaluate characteristics of perturbed cancer-related GRNs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

save_dir <- '../results_new/GRNs'
if(dir.exists(save_dir)){
    unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

fdr <- 0.05

## TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
##----------------------------------------------------------------------

dbd_purt <- data.table::fread('../data/Events_perturbing_DBD.txt')
ed_purt <- data.table::fread('../data/Events_perturbing_ED.txt')
hd_purt <- data.table::fread('../data/Events_perturbing_HD.txt')


##---- Using TFLink ----
tf_ensemb_map <- as.data.frame(data.table::fread('../data/TF_ensembl_uniprot.txt', sep='\t'))
diff_expr_files <- gtools::mixedsort(list.files('../data/Diff_expr', full.names=TRUE))
ensembl_gene_map <- data.table::fread('../data/ensembl_name_map.txt')
##---- grn without direction downloaded from here: https://tflink.net/download/
grn <- as.data.frame(data.table::fread('../data/TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv'))
grn_filt <- grn[grn$`Name.TF` %in% tfs$Gene_Symbol,]
##---- 14,439 edges ----

input_dir <- '../data/PSI_data'
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))
all_files <- all_files[-4]
all_files_raw <- all_files_raw[-4]
all_cancer <- substr(basename(all_files), 1,4)
paired_sam <- data.table::fread('../data/cancer_paired_samples.txt')
paired_sam <- paired_sam[-4]


##--- splicing events affecting TFs ----------------------------------------------------------------------
grns <- list()
grns_uni <- list()
recounts <- c()
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/DBD_perturbed_GRNs.xlsx'))

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    tempdbd1 <- unique(dbd_purt[dbd_purt$CANCER == all_cancer[k], ]$AS)
    tempdbd2 <- unique(ed_purt[ed_purt$CANCER == all_cancer[k], ]$AS)
    tempdbd3 <- unique(hd_purt[hd_purt$CANCER == all_cancer[k], ]$AS)

    tempdbd <- union(union(tempdbd1, tempdbd2), tempdbd3)
    tempy <- temp[temp$as_id %in% tempdbd, ]

    # temptfs <- unique(tempy$symbol)

    # for(j in 1:length(temptfs)){
    #     tempx <- tempy[tempy$symbol == temptfs[j], ]
    #     if(nrow(tempx)>1){
    #         print(j)
    #     }
    # }

    tempy$POSP <- tempy$POS/paired_sam[[2]][k]
    tempy$NEGP <- tempy$NEG/paired_sam[[2]][k]

    temptfs <- unique(tempy$symbol)
    tempgrn1 <- grn_filt[grn_filt$`Name.TF` %in% temptfs, ]
    tempgrn1 <- tempgrn1[,c(1,2,5,6,7,8)]

    ##-- add which event contributes ---
    asid <- c()
    for(j in 1:length(tempgrn1[[1]])){
        tempyu <- tempy[tempy$symbol == tempgrn1$`Name.TF`[j], ]
        tdd <- unlist(mapply(function(x,y) paste0(x,'_',y),tempyu$as_id, tempyu$splice_type))
        # asid <- c(asid, paste(unlist(mapply(function(x,y) paste0(x,'_',y),tdd, signif(tempyu$MEDIAN_DIFF, 3))), collapse=';') )
        asid <- c(asid, paste(unlist(tdd), collapse=';') )
    }
    ##-----------------------------------
    tempgrn1$as_id <- asid
    
    grns[[k]] <- igraph::graph_from_data_frame(tempgrn1[,c(3,4)])
    grns_uni[[k]] <- igraph::graph_from_data_frame(tempgrn1[,c(1,2)])

    recounts <- c(recounts, length(tempgrn1[[1]]))
    ## save excel sheet ----
    tempz <- tempgrn1
    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tempz)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/DBD_perturbed_GRNs.xlsx'), overwrite = T)

}


##--- Number of GR interactions ---
pdata1 <- data.frame(cancer=all_cancer, count=recounts)
p <- ggplot(pdata1, aes(cancer, count)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(pdata1$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of perturbed \ngene regulatory interactions", limits=c(0,maxv+300)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/DBD_ED_perturbed_GRIs.png"),width=3.5, height=3, dpi=400)


##--- Overlap of GR interactions among cancer types ---
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
    temp_unn <- igraph::make_empty_graph(directed=TRUE)
    for(i in 1:loop1){
        temp_ovl <- grns_uni[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- igraph::intersection(temp_ovl, grns_uni[[which(all_cancer == temp_comb[[i]][j])]])
            }
        }
        temp_unn <- igraph::union(temp_unn, temp_ovl)
    }

    num_events_combs[[k]] <- temp_unn
    num_events_combs_counts <- c(num_events_combs_counts, igraph::ecount(temp_unn))
    cat('Cancer',k,'of',length(all_cancer),'done\n')
}

pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)

p <- ggplot(pdata, aes(cancer, count, label=count)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(pdata$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of splicing events", limits=c(0,maxv+300)) +
geom_text(position=position_dodge(0.9), hjust=0, vjust=0, angle=75, size=2.5)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/DBD_ED_perturbed_GRIs_overlap.png"),width=3.5, height=3, dpi=400)


wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/DBD_perturbed_GRNs_atleast.xlsx'))
tempgrnx <- grn_filt[,c(1,2,5,6,7,8)]

for(k in 1:length(all_cancer)){

    tempx <- igraph::as_data_frame(num_events_combs[[k]])

    if(nrow(tempx) != 0){

        ncnt <- c()
        for(j in 1:length(tempx[[1]])){
            wh1 <- which(tempgrnx$`UniprotID.TF` == tempx[[1]][j])
            wh2 <- which(tempgrnx$`UniprotID.Target` == tempx[[2]][j])
            wh <- intersect(wh1, wh2)
            ncnt <- c(ncnt, wh)
        }

        tempz <- tempgrnx[ncnt,]
        openxlsx::addWorksheet(wb1, sheetName = paste0('At least ',k))
        openxlsx::writeData(wb1, sheet = paste0('At least ',k), tempz)
        openxlsx::saveWorkbook(wb1, paste0(save_dir,'/DBD_perturbed_GRNs_atleast.xlsx'), overwrite = T)

    }

}

# ##--- KEGG enrichment of the target genes ----------------------------------------------
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

# temp_grn <- igraph::as_data_frame(num_events_combs[[10]])
# temp_gene <- union(temp_grn[[1]], temp_grn[[2]])
# # num_gene <- length(temp_gene)
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

# #-- save excel file ---
# wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/KEGG_enrichment_10cancer.xlsx'))
# temp1a <- ktdx[order(ktdx$FDR), ]
# colnames(temp1a) <- c('KEGG ID', 'KEGG pathway name', 'pvalue', 'FDR')
# openxlsx::addWorksheet(wb1, sheetName = 'Atleast 10 cancer')
# openxlsx::writeData(wb1, sheet = 'Atleast 10 cancer', temp1a)
# openxlsx::saveWorkbook(wb1, paste0(save_dir,'/KEGG_enrichment_10cancer.xlsx'), overwrite = T)

# ##--------------------------------------------------------------------------
# ##--- save for each cancer type ----
# wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/KEGG_enrichment_per_cancer.xlsx'))

# for(k in 1:length(all_cancer)){

#     pval <- c()
#     Kterm <- c()
#     Ktermn <- c()
#     num_gene <- c()

#     temp_grn <- igraph::as_data_frame(grns_uni[[k]])
#     temp_gene <- union(temp_grn[[1]], temp_grn[[2]])
#     # num_gene <- length(temp_gene)
#     temp_kegg <- allkegg[allkegg$Uniprotswissprot %in% temp_gene, ]
#     sampleSize <- length(unique(temp_kegg$Uniprotswissprot))
#     tempid <- unique(temp_kegg$Go_id)

#     for(j in 1:length(tempid)){
#         setA <- unique(allkegg[allkegg$Go_id == tempid[j], ]$Uniprotswissprot)
#         setB <- unique(temp_kegg[temp_kegg$Go_id == tempid[j], ]$Uniprotswissprot)
#         ## hypergeometric test
#         hyp <- phyper(length(setB)-1,length(setA),bgSize-length(setA), sampleSize,lower.tail = FALSE)
#         pval <- c(pval, hyp)
#         Kterm <- c(Kterm, tempid[j])
#         Ktermn <- c(Ktermn, unique(temp_kegg[temp_kegg$Go_id == tempid[j], ]$Pathways_name))
#     }

#     ktd <- data.frame(Kterm, Ktermn, pval)
#     ktd$FDR <- p.adjust(ktd$pval, 'fdr')
#     ktdx <- ktd[ktd$FDR < 0.05, ]

#     #-- save excel file ---
#     temp1a <- ktdx[order(ktdx$FDR), ]
#     colnames(temp1a) <- c('KEGG ID', 'KEGG pathway name', 'pvalue', 'FDR')
#     openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
#     openxlsx::writeData(wb1, sheet = all_cancer[k], temp1a)
#     openxlsx::saveWorkbook(wb1, paste0(save_dir,'/KEGG_enrichment_per_cancer.xlsx'), overwrite = T)

# }

# ###-------------------------------------------------------------------------------------



# ##--- Correlations between events and gene expressions ---
# store_dir <- '../data/correlations'
# all_sigs <- gtools::mixedsort(list.files(store_dir, pattern='*.txt', full.names=TRUE))
# wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/DBD_perturbed_GRNs_gene_cor.xlsx'))
# wb11 <- openxlsx::createWorkbook(paste0(save_dir,'/DBD_perturbed_GRNs_gene_cor_max.xlsx'))

# recounts <- c()

# for(k in 1:length(all_cancer)){

#     tempcr <- data.table::fread(all_sigs[k])
#     tempcr <- tempcr[tempcr$CFDR < fdr, ]
#     tempcr <- tempcr[tempcr$RFDR < fdr, ]

#     temp <- data.table::fread(all_files[k], sep='\t')
#     temp$POSP <- temp$POS/paired_sam[[2]][k]
#     temp$NEGP <- temp$NEG/paired_sam[[2]][k]

#     tempdbd1 <- unique(dbd_purt[dbd_purt$CANCER == all_cancer[k], ]$AS)
#     tempdbd2 <- unique(ed_purt[ed_purt$CANCER == all_cancer[k], ]$AS)
#     tempdbd <- union(tempdbd1, tempdbd2)
#     tempy <- temp[temp$as_id %in% tempdbd, ]

#     temptfs <- unique(tempy$symbol)
#     tempgrn1 <- grn_filt[grn_filt$`Name.TF` %in% temptfs, ]
#     tempgrn1 <- tempgrn1[,c(1,2,5,6,7,8)]
#     dbd_grn <- igraph::graph_from_data_frame(tempgrn1[,c(3,4)])

#     temptf <- c()
#     posp <- c()
#     negp <- c()
#     mdiff <- c()

#     for(j in 1:length(tempcr[[1]])){
#         temptf <- c(temptf, temp[temp$as_id == tempcr$ASID[j], ]$symbol)
#         posp <- c(posp, temp[temp$as_id == tempcr$ASID[j], ]$POSP)
#         negp <- c(negp, temp[temp$as_id == tempcr$ASID[j], ]$NEGP)
#         mdiff <- c(mdiff, temp[temp$as_id == tempcr$ASID[j], ]$MEDIAN_DIFF)
#     }

#     tempcr$TF <- temptf
#     tempcr$POSP <- posp
#     tempcr$NEGP <- negp
#     tempcr$MEDIAN_DIFF <- mdiff

#     tempint <- tempcr
#     cor <- c()
#     rpval <- c()
#     asid <- c()
#     posp <- c()
#     negp <- c()
#     mdiff <- c()
    
#     for(j in 1:length(tempgrn1[[1]])){
#         wh1 <- which(tempint$TF == tempgrn1$`Name.TF`[j])
#         wh2 <- which(tempint$GENE == tempgrn1$`Name.Target`[j])
#         wh <- intersect(wh1, wh2)
#         if(length(wh) != 0){
#             ct <- c()
#             rt <- c()
#             ast <- c()
#             psp <- c()
#             ngp <- c()
#             miff <- c()
#             for(i in 1:length(wh)){
#                 ct <- c(ct, signif(tempint$COR[wh[i]],3))
#                 rt <- c(rt, signif(tempint$RFDR[wh[i]],3))
#                 ast <- c(ast, tempint$ASID[wh[i]])
#                 psp <- c(psp, signif(tempint$POSP[wh[i]],3))
#                 ngp <- c(ngp, signif(tempint$NEGP[wh[i]],3))
#                 miff <- c(miff, signif(tempint$MEDIAN_DIFF[wh[i]],3))
#             }
#             cor <- c(cor, paste(ct, collapse=';'))
#             rpval <- c(rpval, paste(rt, collapse=';'))
#             asid <- c(asid, paste(ast, collapse=';'))
#             posp <- c(posp, paste(psp, collapse=';'))
#             negp <- c(negp, paste(ngp, collapse=';'))
#             mdiff <- c(mdiff, paste(miff, collapse=';'))
#         }else{
#             cor <- c(cor, NA)
#             rpval <- c(rpval, NA)
#             asid <- c(asid, NA)
#             posp <- c(posp, NA)
#             negp <- c(negp, NA)
#             mdiff <- c(mdiff, NA)
#         }
#     }

#     tempgrn1$COR <- cor
#     tempgrn1$RFDR <- rpval
#     tempgrn1$ASID <- asid
#     tempgrn1$POSP <- posp
#     tempgrn1$NEGP <- negp
#     tempgrn1$MEDIAN_DIFF <- mdiff

#     ## save excel sheet ----
#     tempz <- tempgrn1
#     tempz <- tempz[, -c(5,6)]
#     openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
#     openxlsx::writeData(wb1, sheet = all_cancer[k], tempz)
#     openxlsx::saveWorkbook(wb1, paste0(save_dir,'/DBD_perturbed_GRNs_gene_cor.xlsx'), overwrite = T)

#     tempzz <- tempz[complete.cases(tempz),]
#     recounts <- c(recounts, length(tempzz[[1]]))

#     ##-- max difference ---
#     mdi <- lapply(as.list(tempzz$MEDIAN_DIFF), function(x) as.numeric(unlist(strsplit(x, '[;]'))))
#     maxdiff_pos <- unlist(lapply(mdi, function(x) which(abs(x) == max(abs(x)))[1] ))
#     maxdiff <- unlist(lapply(mdi, function(x) x[which(abs(x) == max(abs(x)))[1]] ))

#     ##-- corresponding cor, rfdr,  asid, posp, negp
#     ncorr <- mapply(function(x,y) as.numeric(unlist(strsplit(x, '[;]'))[y]), as.list(tempzz$COR), as.list(maxdiff_pos))
#     nrfdr <- mapply(function(x,y) as.numeric(unlist(strsplit(x, '[;]'))[y]), as.list(tempzz$RFDR), as.list(maxdiff_pos))
#     nasid <- mapply(function(x,y) unlist(strsplit(x, '[;]'))[y], as.list(tempzz$ASID), as.list(maxdiff_pos))
#     nposp <- mapply(function(x,y) as.numeric(unlist(strsplit(x, '[;]'))[y]), as.list(tempzz$POSP), as.list(maxdiff_pos))
#     nnegp <- mapply(function(x,y) as.numeric(unlist(strsplit(x, '[;]'))[y]), as.list(tempzz$NEGP), as.list(maxdiff_pos))

#     tempzz$NCOR <- ncorr
#     tempzz$NPOSP <- nposp
#     tempzz$NNEGP <- nnegp
#     tempzz$NASID <- nasid
#     tempzz$NRFDR <- nrfdr
#     tempzz$NMEDIAN_DIFF <- maxdiff
#     tempzz <- tempzz[,c(1,2,3,4,11,12,13,14,15,16)]
#     tempzz <- tempzz[order(-abs(tempzz$NMEDIAN_DIFF)),]

#     openxlsx::addWorksheet(wb11, sheetName = all_cancer[k])
#     openxlsx::writeData(wb11, sheet = all_cancer[k], tempzz)
#     openxlsx::saveWorkbook(wb11, paste0(save_dir,'/DBD_perturbed_GRNs_gene_cor_max.xlsx'), overwrite = T)

# }

# ##--- Number of GR interactions with at least one splice event-gene correlations ---
# pdata1$countcor <- recounts
# pdata1$countx <- pdata1$count-pdata1$countcor
# allc <- pdata1$count
# pdata2 <- pdata1[,-2]
# colnames(pdata2) <- c('Cancer','Yes','No')
# pdata <- reshape2::melt(pdata2)
# pdata$count <- rep(allc,2)

# p <- ggplot(pdata, aes(Cancer, value, fill=variable)) + 
# geom_bar(stat="identity",position="stack")+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of GRIs", limits=c(0,max(pdata$count))) +
# # geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# geom_text(aes(label=value), position=position_stack(vjust=0.5), size=2)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill=guide_legend(title="Significant event-\ngene correlations",ncol=1))
# ggsave(p,filename=paste0(save_dir,"/DBD_ED_perturbed_GRIs.png"),width=5, height=3, dpi=400)




