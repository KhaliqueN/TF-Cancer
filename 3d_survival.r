##############################################################################################
# Purpose: Which TF splicing events are survival associated
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

save_dir <- '../results'

## TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
##----------------------------------------------------------------------

input_dir <- '../data/PSI_data'
fdr <- 0.05
diff <- 0.1
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))

all_cancer <- substr(basename(all_files), 1,4)

##-- contradictory event pairs ------------------------------
cont_pairs <- data.table::fread(paste0('../data/contradictory_event_pair_DBD.txt'))

##-- events affecting TFs ------------------------------------
events_tf_pos <- list()
events_tf_neg <- list()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy_pos <- tempy[tempy$MEDIAN_DIFF > 0, ]
    tempy_neg <- tempy[tempy$MEDIAN_DIFF < 0, ]
    events_tf_pos[[k]] <- tempy_pos$as_id
    events_tf_neg[[k]] <- tempy_neg$as_id
}

##--- Overlap with the survival associated events ---------------
## Study: OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers, NAR, 2022
## Downloaded http://47.98.127.64:8080/beta/download?fileName=spliceseq_clinical_as_survival.csv.gz
splicing <- data.table::fread('../data/spliceseq_clinical_as_survival.csv')
splicing <- as.data.frame(splicing)

p1 <- c()
p2 <- c()
md1 <- c()
md2 <- c()
cancer <- c()
events_tf_surv1 <- list()
for(k in 1:length(all_cancer)){

    events_tf <- union(events_tf_pos[[k]], events_tf_neg[[k]])
    temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
    temp1$SE <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
    whi <- intersect(temp1$SE, events_tf)
    tempy <- temp1[temp1$SE %in% whi, ]
    tempy <- tempy[tempy$CI_Type == 'overall_survival',]
    events_tf_surv1[[k]] <- tempy$SE

    ##--- which event pairs show opposite correlations with survival --
    tgenes_c <- plyr::count(tempy$Gene_Symbol)
    tgenes <- tgenes_c$x[which(tgenes_c$freq > 1)]

    if(length(tgenes) != 0){

        for(j in 1:length(tgenes)){
            tempg <- tempy[tempy$Gene_Symbol == tgenes[j], ]
            loop1 <- nrow(tempg)-1
            loop2 <- nrow(tempg)
            for(i in 1:loop1){
                m <- i+1
                for(ii in m:loop2){
                    p1 <- c(p1, paste0(tempg[i,]$Gene_Symbol[1],'_',tempg[i,]$SE,'_',tempg[i,]$Splice_Type))
                    p2 <- c(p2, paste0(tempg[ii,]$Gene_Symbol[1],'_',tempg[ii,]$SE,'_',tempg[ii,]$Splice_Type))
                    md1 <- c(md1, tempg[i, ]$Hazard_Ratio)
                    md2 <- c(md2, tempg[ii, ]$Hazard_Ratio)
                    cancer <- c(cancer, all_cancer[k])
                }
            }
        }
    }

    cat('cancer', k,'of',length(all_cancer),'done\n')

}
tempdata1 <- data.frame(CANCER=cancer, AS1=p1, AS2=p2, MD1=md1, MD2=md2)

p1 <- c()
p2 <- c()
md1 <- c()
md2 <- c()
cancer <- c()
events_tf_surv2 <- list()
for(k in 1:length(all_cancer)){

    events_tf <- union(events_tf_pos[[k]], events_tf_neg[[k]])
    temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
    temp1$SE <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
    whi <- intersect(temp1$SE, events_tf)
    tempy <- temp1[temp1$SE %in% whi, ]
    tempy <- tempy[tempy$CI_Type == 'disease_specific_survival',]
    events_tf_surv2[[k]] <- tempy$SE

    ##--- which event pairs show opposite correlations with survival --
    tgenes_c <- plyr::count(tempy$Gene_Symbol)
    tgenes <- tgenes_c$x[which(tgenes_c$freq > 1)]

    if(length(tgenes) != 0){

        for(j in 1:length(tgenes)){
            tempg <- tempy[tempy$Gene_Symbol == tgenes[j], ]
            loop1 <- nrow(tempg)-1
            loop2 <- nrow(tempg)
            for(i in 1:loop1){
                m <- i+1
                for(ii in m:loop2){
                    p1 <- c(p1, paste0(tempg[i,]$Gene_Symbol[1],'_',tempg[i,]$SE,'_',tempg[i,]$Splice_Type))
                    p2 <- c(p2, paste0(tempg[ii,]$Gene_Symbol[1],'_',tempg[ii,]$SE,'_',tempg[ii,]$Splice_Type))
                    md1 <- c(md1, tempg[i, ]$Hazard_Ratio)
                    md2 <- c(md2, tempg[ii, ]$Hazard_Ratio)
                    cancer <- c(cancer, all_cancer[k])
                }
            }
        }
    }

    cat('cancer', k,'of',length(all_cancer),'done\n')

}
tempdata2 <- unique(data.frame(CANCER=cancer, AS1=p1, AS2=p2, MD1=md1, MD2=md2))

p1 <- c()
p2 <- c()
md1 <- c()
md2 <- c()
cancer <- c()
events_tf_surv3 <- list()
for(k in 1:length(all_cancer)){

    events_tf <- union(events_tf_pos[[k]], events_tf_neg[[k]])
    temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
    temp1$SE <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
    whi <- intersect(temp1$SE, events_tf)
    tempy <- temp1[temp1$SE %in% whi, ]
    tempy <- tempy[tempy$CI_Type == 'disease_free_interval',]
    events_tf_surv3[[k]] <- tempy$SE

    ##--- which event pairs show opposite correlations with survival --
    tgenes_c <- plyr::count(tempy$Gene_Symbol)
    tgenes <- tgenes_c$x[which(tgenes_c$freq > 1)]

    if(length(tgenes) != 0){

        for(j in 1:length(tgenes)){
            tempg <- tempy[tempy$Gene_Symbol == tgenes[j], ]
            loop1 <- nrow(tempg)-1
            loop2 <- nrow(tempg)
            for(i in 1:loop1){
                m <- i+1
                for(ii in m:loop2){
                    p1 <- c(p1, paste0(tempg[i,]$Gene_Symbol[1],'_',tempg[i,]$SE,'_',tempg[i,]$Splice_Type))
                    p2 <- c(p2, paste0(tempg[ii,]$Gene_Symbol[1],'_',tempg[ii,]$SE,'_',tempg[ii,]$Splice_Type))
                    md1 <- c(md1, tempg[i, ]$Hazard_Ratio)
                    md2 <- c(md2, tempg[ii, ]$Hazard_Ratio)
                    cancer <- c(cancer, all_cancer[k])
                }
            }
        }
    }

    cat('cancer', k,'of',length(all_cancer),'done\n')

}
tempdata3 <- unique(data.frame(CANCER=cancer, AS1=p1, AS2=p2, MD1=md1, MD2=md2))

p1 <- c()
p2 <- c()
md1 <- c()
md2 <- c()
cancer <- c()
events_tf_surv4 <- list()
for(k in 1:length(all_cancer)){

    events_tf <- union(events_tf_pos[[k]], events_tf_neg[[k]])
    temp1 <- splicing[splicing$Cancer_Type == all_cancer[k], ]
    temp1$SE <- unlist(lapply(strsplit(temp1$Splice_Event, '[_]'), '[[', 3))
    whi <- intersect(temp1$SE, events_tf)
    tempy <- temp1[temp1$SE %in% whi, ]
    tempy <- tempy[tempy$CI_Type == 'progression_free_interval',]
    events_tf_surv4[[k]] <- tempy$SE

    ##--- which event pairs show opposite correlations with survival --
    tgenes_c <- plyr::count(tempy$Gene_Symbol)
    tgenes <- tgenes_c$x[which(tgenes_c$freq > 1)]

    if(length(tgenes) != 0){
        for(j in 1:length(tgenes)){
            tempg <- tempy[tempy$Gene_Symbol == tgenes[j], ]
            loop1 <- nrow(tempg)-1
            loop2 <- nrow(tempg)
            for(i in 1:loop1){
                m <- i+1
                for(ii in m:loop2){
                    p1 <- c(p1, paste0(tempg[i,]$Gene_Symbol[1],'_',tempg[i,]$SE,'_',tempg[i,]$Splice_Type))
                    p2 <- c(p2, paste0(tempg[ii,]$Gene_Symbol[1],'_',tempg[ii,]$SE,'_',tempg[ii,]$Splice_Type))
                    md1 <- c(md1, tempg[i, ]$Hazard_Ratio)
                    md2 <- c(md2, tempg[ii, ]$Hazard_Ratio)
                    cancer <- c(cancer, all_cancer[k])
                }
            }
        }
    }

    cat('cancer', k,'of',length(all_cancer),'done\n')

}
tempdata4 <- unique(data.frame(CANCER=cancer, AS1=p1, AS2=p2, MD1=md1, MD2=md2))

tempdatax <- rbind(rbind(rbind(tempdata1, tempdata2), tempdata3), tempdata4)
# tempdata <- unique(tempdatax[,c(1,2,3)])

events_tf_surv <- list()
for(k in 1:length(all_cancer)){
    events_tf_surv[[k]] <- union(union(union(events_tf_surv1[[k]], events_tf_surv2[[k]]), events_tf_surv3[[k]]), events_tf_surv4[[k]])
}

##--- Number of events associated with survival ---------------------------------
pdata <- data.frame(cancer=all_cancer, count=lengths(events_tf_surv))
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of events \nassociated with survival", limits=c(0,max(pdata$count)+20)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_survival.png"),width=3.5, height=3, dpi=400)



##-- contradictory event pairs ---
ttdata <- tempdatax
wh1 <- which(ttdata$MD1 < 1 & ttdata$MD2 > 1)
wh2 <- which(ttdata$MD2 < 1 & ttdata$MD1 > 1)
ttdata <- ttdata[union(wh1, wh2), ]
events_tf_surv_opp <- c()
for(k in 1:length(all_cancer)){
    ttgr <- ttdata[ttdata$CANCER == all_cancer[k], ]
    tgr <- igraph::simplify(igraph::graph_from_data_frame(ttgr[,c(2,3)], directed=FALSE))
    events_tf_surv_opp <- c(events_tf_surv_opp, igraph::ecount(tgr))
}

pdata <- data.frame(x=all_cancer, freq=events_tf_surv_opp)#plyr::count(ttdata[[1]])
p <- ggplot(pdata, aes(x, freq)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of contradictory event pairs", limits=c(0,max(pdata$freq)+20)) +
geom_text(aes(label=freq), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0("../results/Survival_contradictory_event_pairs_all.png"),width=3.5, height=3, dpi=400)

data.table::fwrite(ttdata, '../data/contradictory_event_survival.txt', sep='\t', row.names=FALSE, quote=FALSE)

# ##--- Pairs of events affecting the same gene that are contradictory in terms of survival association ----
# wh1 <- which(!is.na(fdata$OS1))
# wh2 <- which(!is.na(fdata$OS2))
# whi <- intersect(wh1, wh2)
# fdatax <- fdata[whi, ]
# wh1 <- which(fdatax$OS1 > 1 & fdatax$OS2 < 1)
# wh2 <- which(fdatax$OS1 < 1 & fdatax$OS2 > 1)
# wh <- union(wh1, wh2)
# fdatay <- fdatax[wh,]

# pdata <- plyr::count(fdatay$CANCER)
# p <- ggplot(pdata, aes(x, freq)) + 
# geom_bar(stat="identity",position=position_dodge())+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of contradictory event pairs \nalso showing contradictory \n survival associations", limits=c(0,120)) +
# geom_text(aes(label=freq), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
# ggsave(p,filename=paste0(save_dir,"/Event_pairs_opp_survival.png"),width=3.5, height=3, dpi=400)


# ##--- overlap of contradictory survival events among cancer types --------------------
# fdatay <- ttdata
# num_sig_cors_d <- list()
# for(k in 1:length(all_cancer)){
#     fdatayt <- fdatay[fdatay$CANCER == all_cancer[k], ]
#     num_sig_cors_d[[k]] <- igraph::graph_from_data_frame(fdatayt[,c(2,3)], directed=FALSE)
# }


# combs <- list() ## store all combinations
# for(k in 1:length(all_cancer)){
#     combs[[k]] <- combn(all_cancer, k)
# }

# num_events_combs <- list()
# num_events_combs_counts <- c()

# for(k in 1:length(all_cancer)){

#     temp_comb <- as.data.frame(combs[[k]])
#     loop1 <- length(temp_comb)
#     loop2 <- length(temp_comb[[1]])
#     temp_unn <- igraph::make_empty_graph(directed=FALSE)
#     for(i in 1:loop1){
#         temp_ovl <- num_sig_cors_d[[which(all_cancer == temp_comb[[i]][1])]]
#         if(loop2 > 1){
#             for(j in 2:loop2){
#                 temp_ovl <- igraph::intersection(temp_ovl, num_sig_cors_d[[which(all_cancer == temp_comb[[i]][j])]])
#             }
#         }
#         temp_unn <- igraph::union(temp_unn, temp_ovl)
#     }

#     num_events_combs[[k]] <- temp_unn
#     num_events_combs_counts <- c(num_events_combs_counts, igraph::ecount(temp_unn))
#     cat('Cancer',k,'of',length(all_cancer),'done\n')
# }

# pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)

# p <- ggplot(pdata, aes(cancer, count, label=count)) + 
# geom_bar(stat="identity")+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# maxv <- max(pdata$count)
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="# of cancer types") + 
# scale_y_continuous(name="# of splicing events", limits=c(0,maxv+100)) +
# geom_text(position=position_dodge(0.9), hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
# ggsave(p,filename=paste0("../results/Contradictory_survival_events_overlap.png"),width=7, height=3, dpi=400)



##------------------------------------------------------------------------------------

##--- Pairs of events affecting the same gene that are contradictory in terms of 
##--- splicing as well as survival association ----
fdatay <- ttdata
num_sig_cors_d <- list()
num_sig_cors_d_count <- c()
for(k in 1:length(all_cancer)){
    fdatayt <- fdatay[fdatay$CANCER == all_cancer[k], ]
    fdatayu <- cont_pairs[cont_pairs$CANCER == all_cancer[k], ]
    fdatayt_gr <- igraph::graph_from_data_frame(fdatayt[,c(2,3)], directed=FALSE)
    fdatayu_gr <- igraph::graph_from_data_frame(fdatayu[,c(2,3)], directed=FALSE)
    num_sig_cors_d[[k]] <- igraph::intersection(fdatayt_gr, fdatayu_gr)
    num_sig_cors_d_count <- c(num_sig_cors_d_count, igraph::ecount(igraph::intersection(fdatayt_gr, fdatayu_gr)))
}


pdata <- data.frame(x=all_cancer, freq=num_sig_cors_d_count)#plyr::count(fdataz$CANCER)
p <- ggplot(pdata, aes(x, freq)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of contradictory event pairs showing \ncontradictory survival associations", 
    limits=c(0,max(pdata$freq)+20)) +
geom_text(aes(label=freq), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Contradictory_event__pairs_opp_survival.png"),width=3.5, height=3, dpi=400)




##--- overlap of contradictory survival events among cancer types --------------------
# fdatay <- ttdata
# num_sig_cors_d <- list()
# for(k in 1:length(all_cancer)){
#     fdatayt <- fdatay[fdatay$CANCER == all_cancer[k], ]
#     num_sig_cors_d[[k]] <- igraph::graph_from_data_frame(fdatayt[,c(2,3)], directed=FALSE)
# }


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
    temp_unn <- igraph::make_empty_graph(directed=FALSE)
    for(i in 1:loop1){
        temp_ovl <- num_sig_cors_d[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- igraph::intersection(temp_ovl, num_sig_cors_d[[which(all_cancer == temp_comb[[i]][j])]])
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
scale_y_continuous(name="# of splicing events", limits=c(0,maxv+20)) +
geom_text(position=position_dodge(0.9), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0("../results/Contradictory_survival_events_overlap.png"),width=7, height=3, dpi=400)
