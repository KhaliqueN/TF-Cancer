##############################################################################################
# Purpose: Evaluate characteristics of perturbed cancer-related GRNs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

save_dir <- '../results_new/GRNs'
if(!dir.exists(save_dir)){
    dir.create(save_dir, recursive=TRUE)
}

fdr <- 0.05

## TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
##----------------------------------------------------------------------

dbd_purt <- data.table::fread('../data/Events_perturbing_DBD.txt')


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

all_cancer <- substr(basename(all_files), 1,4)


##--- splicing events affecting TFs ----------------------------------------------------------------------
grns <- list()
recounts <- c()
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/DBD_perturbed_GRNs.xlsx'))

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    tempdbd <- unique(dbd_purt[dbd_purt$CANCER == all_cancer[k], ]$AS)
    tempy <- temp[temp$as_id %in% tempdbd, ]

    temptfs <- unique(tempy$symbol)
    tempgrn1 <- grn_filt[grn_filt$`Name.TF` %in% temptfs, ]
    tempgrn1 <- tempgrn1[,c(1,2,5,6,7,8)]
    
    grns[[k]] <- igraph::graph_from_data_frame(tempgrn1[,c(3,4)])
    recounts <- c(recounts, length(tempgrn1[[1]]))
    ## save excel sheet ----
    tempz <- tempgrn1
    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tempz)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/DBD_perturbed_GRNs.xlsx'), overwrite = T)

}


##--- Number of GR interactions ---
pdata1 <- data.frame(cancer=all_cancer, count=recounts)
# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity")+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of perturbed \ngene regulatory interactions", limits=c(0,(max(pdata$count))+20)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
# ggsave(p,filename=paste0(save_dir,"/DBD_perturbed_GRIs.png"),width=3.5, height=3, dpi=400)



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
        temp_ovl <- grns[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- igraph::intersection(temp_ovl, grns[[which(all_cancer == temp_comb[[i]][j])]])
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
scale_y_continuous(name="# of splicing events", limits=c(0,maxv+150)) +
geom_text(position=position_dodge(0.9), hjust=0, vjust=0, angle=75, size=2.5)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/DBD_perturbed_GRIs_overlap.png"),width=7, height=3, dpi=400)




# ###--- Number of significant event-gene correlations ----------------------------------------------------------
# all_sigs <- gtools::mixedsort(list.files(store_dir, pattern='*.txt', full.names=TRUE))
# num_sig_cors <- c()
# num_sig_cors_d <- list()
# background_as <- igraph::make_empty_graph(n = 0, directed = FALSE)

# for(k in 1:length(all_cancer)){
#     temp <- data.table::fread(all_sigs[k], sep='\t')
#     temp1 <- temp[temp$RFDR < 0.05, ]
#     temp2 <- temp1[temp1$CFDR < 0.05, ]
#     num_sig_cors <- c(num_sig_cors, length(temp2[[1]]))
#     num_sig_cors_d[[k]] <- igraph::graph_from_data_frame(temp2[,c(1,2)], directed=FALSE)
#     background_as <- igraph::union(igraph::graph_from_data_frame(temp2[,c(1,2)], directed=FALSE), background_as)
# }

# pdata <- data.frame(cancer=all_cancer, count=num_sig_cors)

# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity",position="stack")+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of event-gene pairs\nthat significantly correlate", limits=c(0,max(pdata$count)+250)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')#guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#
# ggsave(p,filename=paste0("../results/Event_gene_pairs_cor.png"),width=3.5, height=3, dpi=400)


# ##-- number of event-gene correlations where the event perturbs a DBD ---
# num_sig_corsx <- c()
# num_sig_cors_dx <- list()
# background_asx <- igraph::make_empty_graph(n = 0, directed = FALSE)

# for(k in 1:length(all_cancer)){
#     temp <- data.table::fread(all_sigs[k], sep='\t')
#     temp1 <- temp[temp$RFDR < 0.05, ]
#     temp2 <- temp1[temp1$CFDR < 0.05, ]
#     tempdbd <- unique(dbd_purt[dbd_purt$CANCER == all_cancer[k], ]$AS)
#     temp3 <- temp2[temp2$ASID %in% tempdbd, ]
#     num_sig_corsx <- c(num_sig_corsx, length(temp3[[1]]))
#     num_sig_cors_dx[[k]] <- igraph::graph_from_data_frame(temp3[,c(1,2)], directed=FALSE)
#     background_asx <- igraph::union(igraph::graph_from_data_frame(temp3[,c(1,2)], directed=FALSE), background_as)
# }

# pdata <- data.frame(cancer=all_cancer, count=num_sig_cors)



##--- Correlations between events and gene expressions ---
store_dir <- '../data/correlations'
all_sigs <- gtools::mixedsort(list.files(store_dir, pattern='*.txt', full.names=TRUE))
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/DBD_perturbed_GRNs_gene_cor.xlsx'))
recounts <- c()

for(k in 1:length(all_cancer)){

    tempcr <- data.table::fread(all_sigs[k])
    tempcr <- tempcr[tempcr$CFDR < fdr, ]
    tempcr <- tempcr[tempcr$RFDR < fdr, ]

    temp <- data.table::fread(all_files[k], sep='\t')
    tempdbd <- unique(dbd_purt[dbd_purt$CANCER == all_cancer[k], ]$AS)
    tempy <- temp[temp$as_id %in% tempdbd, ]
    # tempint <- tempcr[tempcr$ASID %in% tempy$as_id, ]
    
    # tempyy <- tempy[tempy$as_id %in% tempint$ASID, ]
    temptfs <- unique(tempy$symbol)
    tempgrn1 <- grn_filt[grn_filt$`Name.TF` %in% temptfs, ]
    tempgrn1 <- tempgrn1[,c(1,2,5,6,7,8)]
    dbd_grn <- igraph::graph_from_data_frame(tempgrn1[,c(3,4)])

    # temptf <- c()
    # for(j in 1:length(tempint[[1]])){
    #     temptf <- c(temptf, tempy[tempy$as_id == tempint$ASID[j], ]$symbol)
    # }
    # tempint$TF <- temptf

    cor <- c()
    rpval <- c()
    asid <- c()
    temptf <- c()
    for(j in 1:length(tempcr[[1]])){
        temptf <- c(temptf, temp[temp$as_id == tempcr$ASID[j], ]$symbol)
    }
    tempcr$TF <- temptf
    tempint <- tempcr
    
    for(j in 1:length(tempgrn1[[1]])){
        wh1 <- which(tempint$TF == tempgrn1$`Name.TF`[j])
        wh2 <- which(tempint$GENE == tempgrn1$`Name.Target`[j])
        wh <- intersect(wh1, wh2)
        if(length(wh) != 0){
            ct <- c()
            rt <- c()
            ast <- c()
            for(i in 1:length(wh)){
                ct <- c(ct, signif(tempint$COR[wh[i]],3))
                rt <- c(rt, signif(tempint$RFDR[wh[i]],3))
                ast <- c(ast, tempint$ASID[wh[i]])
            }
            cor <- c(cor, paste(ct, collapse=';'))
            rpval <- c(rpval, paste(rt, collapse=';'))
            asid <- c(asid, paste(ast, collapse=';'))
        }else{
            cor <- c(cor, NA)
            rpval <- c(rpval, NA)
            asid <- c(asid, NA)
        }
    }

    tempgrn1$COR <- cor
    tempgrn1$RFDR <- rpval
    tempgrn1$ASID <- asid

    ## save excel sheet ----
    tempz <- tempgrn1
    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tempz)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/DBD_perturbed_GRNs_gene_cor.xlsx'), overwrite = T)

    wh <- length(which(!is.na(tempgrn1$ASID)))
    recounts <- c(recounts, wh)
    # cor_grn <- igraph::graph_from_data_frame(tempint[,c(8,2)])
    # bt_grn <- igraph::intersection(cor_grn, dbd_grn)
    # print(igraph::ecount(bt_grn))
    # print(igraph::ecount(bt_grn)/igraph::ecount(cor_grn))

}

##--- Number of GR interactions with at least one splice event-gene correlations ---
pdata1$countcor <- recounts
# colnames(pdata1) <- c('cancer','count','countcor')
# pdata1 <- pdata1[,c(1,2,3,6)]
pdata1$countx <- pdata1$count-pdata1$countcor
allc <- pdata1$count
pdata2 <- pdata1[,-2]
colnames(pdata2) <- c('Cancer','Yes','No')
pdata <- reshape2::melt(pdata2)
pdata$count <- rep(allc,2)

p <- ggplot(pdata, aes(Cancer, value, fill=variable)) + 
geom_bar(stat="identity",position="stack")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of GRIs", limits=c(0,max(pdata$count))) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
geom_text(aes(label=value), position=position_stack(vjust=0.5), size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Significant event-\ngene correlations",ncol=1))
ggsave(p,filename=paste0(save_dir,"/DBD_perturbed_GRIs.png"),width=5, height=3, dpi=400)

# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity")+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of perturbed \ngene regulatory interactions", limits=c(0,(max(pdata$count))+10)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
# ggsave(p,filename=paste0(save_dir,"/DBD_perturbed_GRIs_sig_event-gene_corr.png"),width=3.5, height=3, dpi=400)








# ##--- Motif characterization of cancer related GRNs ---
# ##-- run mfinder on LUSC ------
# ## mfinder downloaded from https://www.weizmann.ac.il/mcb/alon/download/network-motif-software
# ## main.c corrected --> added "extern" in the Role_members variable
# ## Once compiled, add path to ~/.profile
# save_dir <- '../data/motifs'
# if(!dir.exists(save_dir)){
#     dir.create(save_dir)
# }

# for(k in 1:length(all_cancer)){

#     temp <- igraph::as_data_frame(grns[[k]])
#     temp$wt <- rep(1, length(temp[[1]]))
#     data.table::fwrite(temp,'../data/grn.tmp',sep=',', col.names=FALSE, row.names=FALSE, quote=FALSE)
#     cmd <- paste('mfinder ../data/grn.tmp -s 3 -r 10 -f', paste0(save_dir,'/',all_cancer[k]))
#     system(cmd)

# }






