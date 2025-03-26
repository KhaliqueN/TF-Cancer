##############################################################################################
# Purpose: For each cancer type, see whether DBDs are removed by AS
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)
library(dorothea)

save_dir <- '../results'
psi_input <- '../data/PSI_data'
as_input <- '../data/uniprot_Ensembl_Exon_map_DBD_AS'

fdr <- 0.05
diff <- 0.1
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
tf_ensemb_map <- as.data.frame(data.table::fread('../data/TF_ensembl_uniprot.txt', sep='\t'))
tcga_map <- data.table::fread(paste0(psi_input,'/TCGA_SpliceSeq_Gene_Structure.txt'))

all_files <- gtools::mixedsort(list.files(psi_input, pattern='*filtered_PSI.txt', full.names=TRUE))
all_cancer <- substr(basename(all_files), 1,4)
## Some of the event coordinates from TCGA splice seq does not overlap to any of the exons mapped from the canonical protein
## sequence from Uniprot. This is because the canonical seqeunce is not always the longest.



# ##--- Number of TFs with DBDs affected by ES/AP splicing events -------------------------------
# tfs_curated <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
# ## 1340 out of the 1639 curated TFs have at least one known DBD information ---

# ##--- store TFs wth DBD info from uniprot ---
# TFs_with_DBD <- c()
# all_filesx <- list.files(input_dirx, pattern='*.txt', full.names=TRUE)

# for(k in 1:length(all_filesx)){

#     temp <- data.table::fread(all_filesx[k])
#     wh <- which(temp$DBD != '-')
#     if(length(wh) != 0){

#         tid <- strsplit(basename(all_filesx[k]), '[.]')[[1]][1]
#         eid <- tf_ensemb_map[which(tf_ensemb_map$Uniprotswissprot == tid),]$Ensembl_gene_id
#         TFs_with_DBD <- union(TFs_with_DBD, tfs_curated[which(tfs$Ensembl_Gene_ID == eid), ]$Gene_Symbol)
#     }

# }


# ##------------------------------------------
# num_tfs <- c()
# num_tf_events <- c()
# all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI.txt', full.names=TRUE))

# for(k in 1:length(all_cancer)){

#     temp <- data.table::fread(all_files[k], sep='\t')
#     wh1 <- which(temp$FDR < fdr)
#     wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
#     wh <- intersect(wh1, wh2)
#     tempx <- temp[wh, ]
#     tempy <- tempx[which(toupper(tempx$symbol) %in% tfs$Gene_Symbol), ]
#     tempz <- tempy[tempy$splice_type %in% c('ES','AP'), ]
#     all_tfs <- unique(tempz$symbol)

#     ##-- which of the Tfs have at least one DBD info in the uniprot ----
#     cancer_tfs <- intersect(all_tfs, TFs_with_DBD)

#     num_tfs <- c(num_tfs, length(cancer_tfs))
#     num_tf_events <- c(num_tf_events, length(tempz$as_id))
    
# }


# ##--- plot the number of TFs -----------------
# pdata <- data.frame(cancer=all_cancer, count=num_tfs)
# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity",position=position_dodge())+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of TFs affected by a \n significant splicing event with DBDs", limits=c(0,150)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
# ggsave(p,filename=paste0(save_dir,"/Affected_TFs_with_DBDs.png"),width=3.5, height=3, dpi=400)

# ##--- plot the number of splicing events affecting TFs -----------------
# pdata <- data.frame(cancer=all_cancer, count=num_tf_events)
# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity",position=position_dodge())+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of significant events affecting a \n TF with DBDs", limits=c(0,300)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
# ggsave(p,filename=paste0(save_dir,"/Events_affecting_TFs_with_DBDs.png"),width=3.5, height=3, dpi=400)


##----------- Plots -------

ES_affecting_TFs_with_DBD <- vector(mode = "list", length = length(all_cancer))
AP_affecting_TFs_with_DBD <- vector(mode = "list", length = length(all_cancer))
AT_affecting_TFs_with_DBD <- vector(mode = "list", length = length(all_cancer))
AD_affecting_TFs_with_DBD <- vector(mode = "list", length = length(all_cancer))
AA_affecting_TFs_with_DBD <- vector(mode = "list", length = length(all_cancer))
ME_affecting_TFs_with_DBD <- vector(mode = "list", length = length(all_cancer))

TFs_with_affected_DBD_ES <- vector(mode = "list", length = length(all_cancer))
TFs_with_affected_DBD_AP <- vector(mode = "list", length = length(all_cancer))
TFs_with_affected_DBD_AT <- vector(mode = "list", length = length(all_cancer))
TFs_with_affected_DBD_AD <- vector(mode = "list", length = length(all_cancer))
TFs_with_affected_DBD_AA <- vector(mode = "list", length = length(all_cancer))
TFs_with_affected_DBD_ME <- vector(mode = "list", length = length(all_cancer))

for(k in 1:length(all_cancer)){

    all_files <- list.files(as_input, pattern=paste0('*',all_cancer[k],'.txt'), full.names=TRUE)

    all_filesx1 <- c()
    all_filesx2 <- c()
    all_filesx3 <- c()
    all_filesx4 <- c()
    all_filesx5 <- c()
    all_filesx6 <- c()

    events1 <- c()
    events2 <- c()
    events3 <- c()
    events4 <- c()
    events5 <- c()
    events6 <- c()

    for(j in 1:length(all_files)){

        temp <- data.table::fread(all_files[j])
        temp1 <- temp[temp$ES != '-', ]
        temp2 <- temp[temp$AP != '-', ]
        temp3 <- temp[temp$AT != '-', ]
        temp4 <- temp[temp$AD != '-', ]
        temp5 <- temp[temp$AA != '-', ]
        temp6 <- temp[temp$ME != '-', ]

        if(nrow(temp1) > 0){
            if(length(which(temp1$DBD != '-')) > 0){
                all_filesx1 <- c(all_filesx1, all_files[j])
                events1 <- c(events1, setdiff(unique(temp1$ES), '-'))
            }
        }

        if(nrow(temp2) > 0){
            if(length(which(temp2$DBD != '-')) > 0){
                all_filesx2 <- c(all_filesx2, all_files[j])
                events2 <- c(events2, setdiff(unique(temp2$AP), '-'))
            }
        }

        if(nrow(temp3) > 0){
            if(length(which(temp3$DBD != '-')) > 0){
                all_filesx3 <- c(all_filesx3, all_files[j])
                events3 <- c(events3, setdiff(unique(temp3$AT), '-'))
            }
        }

        if(nrow(temp4) > 0){
            if(length(which(temp4$DBD != '-')) > 0){
                all_filesx4 <- c(all_filesx4, all_files[j])
                events4 <- c(events4, setdiff(unique(temp4$AD), '-'))
            }
        }

        if(nrow(temp5) > 0){
            if(length(which(temp5$DBD != '-')) > 0){
                all_filesx5 <- c(all_filesx5, all_files[j])
                events5 <- c(events5, setdiff(unique(temp5$AA), '-'))
            }
        }

        if(nrow(temp6) > 0){
            if(length(which(temp6$DBD != '-')) > 0){
                all_filesx6 <- c(all_filesx6, all_files[j])
                events6 <- c(events6, setdiff(unique(temp6$ME), '-'))
            }
        }

    }

    if(length(all_filesx1) != 0){
        TFs_with_affected_DBD_ES[[k]] <- unlist(lapply(strsplit(basename(all_filesx1), '[_]'),'[[',1))
    }
    if(length(all_filesx2) != 0){
        TFs_with_affected_DBD_AP[[k]] <- unlist(lapply(strsplit(basename(all_filesx2), '[_]'),'[[',1)) 
    }
    if(length(all_filesx3) != 0){
        TFs_with_affected_DBD_AT[[k]] <- unlist(lapply(strsplit(basename(all_filesx3), '[_]'),'[[',1))  
    }
    if(length(all_filesx4) != 0){
        TFs_with_affected_DBD_AD[[k]] <- unlist(lapply(strsplit(basename(all_filesx4), '[_]'),'[[',1))  
    }
    if(length(all_filesx5) != 0){
        TFs_with_affected_DBD_AA[[k]] <- unlist(lapply(strsplit(basename(all_filesx5), '[_]'),'[[',1))  
    }
    if(length(all_filesx6) != 0){
        TFs_with_affected_DBD_ME[[k]] <- unlist(lapply(strsplit(basename(all_filesx6), '[_]'),'[[',1))  
    }

    if(!is.null(events1)){ES_affecting_TFs_with_DBD[[k]] <- events1}
    if(!is.null(events2)){AP_affecting_TFs_with_DBD[[k]] <- events2}
    if(!is.null(events3)){AT_affecting_TFs_with_DBD[[k]] <- events3}
    if(!is.null(events4)){AD_affecting_TFs_with_DBD[[k]] <- events4}
    if(!is.null(events5)){AA_affecting_TFs_with_DBD[[k]] <- events5}
    if(!is.null(events6)){ME_affecting_TFs_with_DBD[[k]] <- events6}

}

TFs_list <- list(TFs_with_affected_DBD_ES, TFs_with_affected_DBD_AP, TFs_with_affected_DBD_AT,
    TFs_with_affected_DBD_AD, TFs_with_affected_DBD_AA, TFs_with_affected_DBD_ME)

Events_list <- list(ES_affecting_TFs_with_DBD, AP_affecting_TFs_with_DBD, AT_affecting_TFs_with_DBD,
    AD_affecting_TFs_with_DBD, AA_affecting_TFs_with_DBD, ME_affecting_TFs_with_DBD)


#----- All TFs and all events ---
all_tfs <- list()
all_events <- list()
for(k in 1:length(all_cancer)){
    temptf <- c()
    tempev <- c()
    for(j in 1:length(TFs_list)){
        temptf <- union(temptf, TFs_list[[j]][[k]])
        tempev <- union(tempev, Events_list[[j]][[k]])
    }
    all_tfs[[k]] <- temptf
    all_events[[k]] <- tempev
}

##--- Events plots ------
Event_count <- c()
as_tags <- c('ES','AP','AT','AD','AA','ME')
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()
for(k in 1:length(all_cancer)){

    tag <- c()
    ect <- c()
    for(j in 1:length(Events_list)){
        ect <- union(ect, Events_list[[j]][[k]])
        tag <- c(tag, rep(as_tags[j], length(Events_list[[j]][[k]])))
    }
    Event_count <- c(Event_count, length(ect))
    temp_count <- plyr::count(tag)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
    
}
##--- plot the number of splicing events affecting TFs -----------------
pdata <- data.frame(cancer=all_cancer, count=Event_count)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of AS events perturbing \na DBD region of TFs", limits=c(0,75)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Events_perturbing_DBDs.png"),width=3.5, height=3, dpi=400)


##--- plot the number of splicing events of different types affecting TFs -----------------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, X=Event_count, Cancer=all_cancer)

pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of significant AS events \naffecting DBDs of TFs")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
guides(fill=guide_legend(title="Alternative splicing type"))
ggsave(ppx,filename=paste0(save_dir, "/Events_types_perturbing_DBDs.png"),width=7, height=3, dpi=400)


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
    tempy <- data.frame(matrix(ncol=0, nrow=0))

    for(i in 1:loop1){
        temp_ovl <- all_events[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- intersect(temp_ovl, all_events[[which(all_cancer == temp_comb[[i]][j])]])
            }
        }
        temp_unn <- union(temp_unn, temp_ovl)
    }
    num_events_combs[[k]] <- temp_unn
    num_events_combs_counts <- c(num_events_combs_counts, length(temp_unn))
}

pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of splicing events", limits=c(0,(max(pdata$count))+10)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Events_perturbing_DBDs_overlap.png"),width=7, height=3, dpi=400)


# ##--- Overlap with the survival associated events ---------------
# ## Study: OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers, NAR, 2022
# ## Downloaded http://47.98.127.64:8080/beta/download?fileName=spliceseq_clinical_as_survival.csv.gz
# splicing <- data.table::fread('../data/spliceseq_clinical_as_survival.csv')
# splicing <- as.data.frame(splicing)

# for(k in 1:length(all_cancer)){
#     temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
#     temp1$SE <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
#     whi <- intersect(temp1$SE, all_events[[k]])
#     tempy <- temp1[temp1$SE %in% whi, ]
 
# }


##--- plot the number of TFs -----------------
pdata <- data.frame(cancer=all_cancer, count=lengths(all_tfs))
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of TFs with DBDs \naffected by a splicing event", limits=c(0,75)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/TFs_with_affected_DBDs.png"),width=3.5, height=3, dpi=400)



##--- which DBD types are affected --------------------------------------
affected_dbds <- list()
for(k in 1:length(all_cancer)){

    dbds <- c()
    for(j in 1:length(all_tfs[[k]])){
        temp <- data.table::fread(paste0(as_input,'/',all_tfs[[k]][j],'_',all_cancer[k],'.txt'))
        temp1 <- temp[temp$DBD != '-', ]
        es <- which(temp1$ES != '-')
        ap <- which(temp1$AP != '-')
        at <- which(temp1$AT != '-')
        ad <- which(temp1$AD != '-')
        aa <- which(temp1$AA != '-')
        me <- which(temp1$ME != '-')
        trem <- union(union(union(union(union(es,ap),at),ad),aa),me)
        dbds <- c(dbds,unique(temp1[trem,]$DBD))
    }

    affected_dbds[[k]] <- plyr::count(dbds)
}

##--- Which DBD types are never affected? Why? --------------------------








#####----- GRNs --------------------------------------------------------

##---- Using TFLink ----
diff_expr_files <- gtools::mixedsort(list.files('../data/Diff_expr', full.names=TRUE))
ensembl_gene_map <- data.table::fread('../data/ensembl_name_map.txt')

##--- grn without direction downloaded from here: https://tflink.net/download/
grn <- as.data.frame(data.table::fread('../data/TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv'))
grn_filt <- grn[grn$`Name.TF` %in% tfs$Gene_Symbol,]
##--- 14,439 edges ----

# ##-- KEGG pathways analysis of the affected TFs -----
# ##--- store TFs wth DBD info from uniprot ---
# TFs_with_DBD <- c()
# all_filesx <- list.files(input_dirx, pattern='*.txt', full.names=TRUE)

# for(k in 1:length(all_filesx)){

#     temp <- data.table::fread(all_filesx[k])
#     wh <- which(temp$DBD != '-')
#     if(length(wh) != 0){
#         tid <- strsplit(basename(all_filesx[k]), '[.]')[[1]][1]
#         TFs_with_DBD <- union(TFs_with_DBD, tid)
#     }
# }

# all_kegg <- data.table::fread('../data/all_human_pathways.txt', sep='\t')
# colnames(all_kegg) <- c('Ensembl_gene_id','Uniprotswissprot','Entrez_id','Go_id','Pathways_name')
# all_kegg1 <- all_kegg[all_kegg$Uniprotswissprot %in% TFs_with_DBD, ]
# all_path <- plyr::count(all_kegg1$Go_id) ## 347 total KEGG pathways
# wh_path <- unique(all_path[which(all_path$freq > 2), ]$x)
# allkeggx <- all_kegg1[all_kegg1$Go_id %in% wh_path, ] ### all pathways remain
# bgSize <- length(unique(allkeggx$Uniprotswissprot))
# pval <- c()
# tcancer <- c()
# Kterm <- c()
# Ktermn <- c()
# num_gene <- c()

# for(k in 1:length(all_cancer)){

#     # ##---- expression matrix ----
#     # temp_expr <- data.table::fread(diff_expr_files[k])
#     # temp_expr$Ensembl_gene_id <- unlist(lapply(strsplit(temp_expr$Ensembl_gene_id, '[.]'),'[[',1))

#     temp_kegg <- allkeggx[allkeggx$Uniprotswissprot %in% all_tfs[[k]], ]
#     sampleSize <- length(unique(temp_kegg$Uniprotswissprot))
#     tempid <- unique(temp_kegg$Go_id)

#     for(j in 1:length(tempid)){
#         setA <- unique(allkeggx[allkeggx$Go_id == tempid[j], ]$Uniprotswissprot)
#         setB <- unique(temp_kegg[temp_kegg$Go_id == tempid[j], ]$Uniprotswissprot)
#         ## hypergeometric test
#         hyp <- phyper(length(setB)-1,length(setA),bgSize-length(setA), sampleSize,lower.tail = FALSE)
#         pval <- c(pval, hyp)
#         Kterm <- c(Kterm, tempid[j])
#         Ktermn <- c(Ktermn, unique(temp_kegg[temp_kegg$Go_id == tempid[j], ]$Pathways_name))
#         tcancer <- c(tcancer, all_cancer[k])
#     }

#     cat('Cancer',k,'of',length(all_cancer),'done\n')
# }

# ktd <- data.frame(tcancer, Kterm, Ktermn, pval)
# ktd$FDR <- p.adjust(ktd$pval, 'fdr')
# ktdx <- ktd[ktd$FDR < 0.05, ] --> nothing enriched


##-- KEGG pathways analysis of the affected genes -----
all_kegg <- data.table::fread('../data/all_human_pathways.txt', sep='\t')
colnames(all_kegg) <- c('Ensembl_gene_id','Uniprotswissprot','Entrez_id','Go_id','Pathways_name')
all_path <- plyr::count(all_kegg$Go_id) ## 347 total KEGG pathways
wh_path <- unique(all_path[which(all_path$freq > 2), ]$x)
allkegg <- all_kegg[all_kegg$Go_id %in% wh_path, ] ### all pathways remain
bgSize <- length(unique(allkegg$Uniprotswissprot))
pval <- c()
tcancer <- c()
Kterm <- c()
Ktermn <- c()
num_gene <- c()

for(k in 1:length(all_cancer)){

    ##---- expression matrix ----
    temp_expr <- data.table::fread(diff_expr_files[k])
    temp_expr$Ensembl_gene_id <- unlist(lapply(strsplit(temp_expr$Ensembl_gene_id, '[.]'),'[[',1))

    temp_gene <- grn_filt[grn_filt$`UniprotID.TF` %in% all_tfs[[k]], ]$`UniprotID.Target`
    # temp_gene_en <- unique(ensembl_gene_map[ensembl_gene_map$Uniprotswissprot %in% temp_gene, ]$Ensembl_gene_id)
    num_gene <- c(num_gene, length(temp_gene))

    temp_kegg <- allkegg[allkegg$Uniprotswissprot %in% temp_gene, ]
    sampleSize <- length(unique(temp_kegg$Uniprotswissprot))
    tempid <- unique(temp_kegg$Go_id)

    for(j in 1:length(tempid)){
        setA <- unique(allkegg[allkegg$Go_id == tempid[j], ]$Uniprotswissprot)
        setB <- unique(temp_kegg[temp_kegg$Go_id == tempid[j], ]$Uniprotswissprot)
        ## hypergeometric test
        hyp <- phyper(length(setB)-1,length(setA),bgSize-length(setA), sampleSize,lower.tail = FALSE)
        pval <- c(pval, hyp)
        Kterm <- c(Kterm, tempid[j])
        Ktermn <- c(Ktermn, unique(temp_kegg[temp_kegg$Go_id == tempid[j], ]$Pathways_name))
        tcancer <- c(tcancer, all_cancer[k])
    }

    cat('Cancer',k,'of',length(all_cancer),'done\n')
}

ktd <- data.frame(tcancer, Kterm, Ktermn, pval)
ktd$FDR <- p.adjust(ktd$pval, 'fdr')
ktdx <- ktd[ktd$FDR < 0.05, ]

#-- save excel file ---
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/KEGG_enrichment.xlsx'))
allKE <- ktdx
for(k in 1:length(all_cancer)){

    temp <- allKE[allKE$tcancer == all_cancer[k], ]
    temp1a <- temp[order(temp$FDR), ]

    colnames(temp1a) <- c('Cancer','KEGG ID', 'KEGG pathway name', 'pvalue', 'FDR')
    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], temp1a)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/KEGG_enrichment.xlsx'), overwrite = T)

}


####----- Pairwise overlaps of enriched KEGG terms ---------------------------------------------------------------------------------
path_bgSize <- length(unique(allkegg$Go_id))
cancer_type <- all_cancer
c1 <- c()
c2 <- c()
value <- c()
pvalue <- c()
pvalueu <- c()
loop1 <- length(cancer_type)-1
loop2 <- length(cancer_type)

for(k in 1:loop1){

    tempx1 <- allKE[allKE$tcancer == cancer_type[k], ]
    m  <- k+1
    for(j in m:loop2){

        # if(k != j){
        tempx2 <- allKE[allKE$tcancer == cancer_type[j], ]
        e1 <- unique(tempx1$Kterm)
        e2 <- unique(tempx2$Kterm)

        ovr <- intersect(e1, e2)
        ur <- union(e1, e2)
        pv <- phyper(length(ovr)-1, length(e1), path_bgSize-length(e1), length(e2), lower.tail=FALSE)
        pvu <- phyper(length(ovr), length(e1), path_bgSize-length(e1), length(e2), lower.tail=TRUE)

        if(length(e1) < length(e2)){
            tempo <- length(ovr)/length(e1)
        }else{
            tempo <- length(ovr)/length(e2)
        }
        value <- c(value, length(ovr))

        pvalue <- c(pvalue, pv)
        pvalueu <- c(pvalueu, pvu)
        c1 <- c(c1, paste0(cancer_type[k],' (',length(e1),')'))
        c2 <- c(c2, paste0(cancer_type[j],' (',length(e2),')'))
    # }
    }
}

qval <- signif(p.adjust(pvalue,'fdr'),2)
qvalu <- signif(p.adjust(pvalueu,'fdr'),2)
qvalx <- rep('Overlap as \nexpected by chance',length(pvalue))
for(i in 1:length(pvalue)){
    if(qval[i] < 0.05){
        qvalx[i] <- 'Higher overlap than\nexpected by chance'
    }
    if(qvalu[i] < 0.05){
        qvalx[i] <- 'Lower overlap than\nexpected by chance'
    }
}

pdata <- data.frame(c1=c1, c2=c2, val=value, pval=qvalx)

cols <- rev(c('#377eb8','#e41a1c'))

p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pval),colour = "white")+
  theme(legend.text=element_text(size=8))+scale_fill_manual(values=cols,drop=FALSE)
basesize <- 8
p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type", y = "Cancer type") +
  scale_y_discrete() +
  scale_x_discrete()+coord_fixed()+
  guides(fill=guide_legend(title="Category"))+
  geom_text(data=pdata,aes(y=c2,label=val),size=2.5)+
  theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0.5, colour = "black"))#+
  # guides(fill='none')
ggsave(p,filename=paste0(save_dir,'/KEGG_overlap.png'),width=4.5, height=3, dpi=300)


###---- unique KEGG terms -------
unique_KE <- list()
uKEs <- c()
for(k in 1:length(cancer_type)){
    counter <- k
    temp_uni <- unique(allKE[allKE$tcancer == cancer_type[k], ]$Ktermn)

    for(j in 1:length(cancer_type)){
        if(counter != j){
            temp_uni <- setdiff(temp_uni, allKE[allKE$tcancer == cancer_type[j], ]$Ktermn)
        }
    }
    unique_KE[[k]] <- temp_uni
    uKEs <- c(uKEs, length(temp_uni))
}

tempg <- data.frame(cancer=cancer_type, count=uKEs)
p <- ggplot(tempg, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of uniquely enriched \nKEGG pathways", limits=c(0,(max(tempg$count))+2)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0.5, vjust=0, angle=0, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Unique_KEGG.png"),width=3.5, height=3, dpi=400)



# #####----- Using Dorothea ----
# all_map <- data.table::fread('../data/ensembl_name_map.txt')
# conf_grn <- as.data.frame(dorothea_hs[dorothea_hs$confidence == 'A', ])
# ##------------- 5,434 edges are confident
# temp_tf <- c()
# temp_target <- c()
# temp_tf_en <- c()
# temp_target_en <- c()
# mult <- c()
# for(k in 1:length(conf_grn[[1]])){

#     ttf <- setdiff(unique(all_map$Uniprotswissprot[which(all_map$HGNC_symbol == conf_grn$tf[k])]),'')
#     if(length(ttf) > 1){
#         mult <- c(mult, k)
#     }
#     temp_tf <- c(temp_tf, ttf[1])

#     ttfx <- setdiff(unique(all_map$Uniprotswissprot[which(all_map$HGNC_symbol == conf_grn$target[k])]),'')
#     if(length(ttfx) > 1){
#         mult <- c(mult, k)
#     }
#     temp_target <- c(temp_target, ttfx[1])

#     ttf <- setdiff(unique(all_map$Ensembl_gene_id[which(all_map$HGNC_symbol == conf_grn$tf[k])]),'')
#     if(length(ttf) > 1){
#         mult <- c(mult, k)
#     }
#     temp_tf_en <- c(temp_tf_en, ttf[1])

#     ttfx <- setdiff(unique(all_map$Ensembl_gene_id[which(all_map$HGNC_symbol == conf_grn$target[k])]),'')
#     if(length(ttfx) > 1){
#         mult <- c(mult, k)
#     }
#     temp_target_en <- c(temp_target_en, ttfx[1])

#     cat('Protein',k, 'of',length(conf_grn[[1]]), 'done\n')
# }


# conf_grn$TF_uniprot <- temp_tf
# conf_grn$Target_uniprot <- temp_target
# conf_grn$TF_Ensembl <- temp_tf_en
# conf_grn$Target_Ensembl <- temp_target_en
# conf_grn1 <- conf_grn[complete.cases(conf_grn), ]
# ###---------- 5,392 edges remain after identifier mappings

# ## number of negative/positive regutions
# neg_num <- which(conf_grn1$mor == -1)
# pos_num <- which(conf_grn1$mor == 1)

# ##-- number of neg / pos regulation disruptions in union of cancer types --
# alltfs <- c()
# for(k in 1:length(all_cancer)){
#     alltfs <- union(alltfs, all_tfs[[k]])
# }
# ## 115 TFs perturbed in at least one cancer type

# temp_grn <- conf_grn1[conf_grn1$TF_uniprot %in% alltfs, ]
# ## Only two of them are present in the confident GRN covering 52 of the 5,392 edges ----









