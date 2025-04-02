##############################################################################################
# Purpose: Plots to see 
# - The number of GENES with splicing events significantly different between normal and cancer tissues
# - Number of TFs with AS events 
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)

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

num_sig_events <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    tempx1 <- unique(tempx$symbol)
    num_sig_events <- c(num_sig_events, length(tempx1))

}

##--- plot the number of splicing events -----------------
pdata <- data.frame(cancer=all_cancer, count=num_sig_events)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of genes affected by \nperturbed splicing events", limits=c(0,(max(pdata$count))+1000)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Sig_events_genes.png"),width=3.5, height=3, dpi=400)


##--- splicing events affecting TFs ----------------------------------------------------------------------
num_tf_events <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    num_tf_events <- c(num_tf_events, length(unique(tempy$symbol)))
    
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
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_genes.png"),width=3.5, height=3, dpi=400)


##-------- Pairwise overlap of slicing events affecting TFs across cancer types --------------------------------------------
TF_splicing_events <- list()
background_as <- c()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    whx <- which(tempx$symbol %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    TF_splicing_events[[k]] <- unique(tempy$symbol)
    background_as <- union(background_as, unique(tempy$symbol))
}

##--- figure for the overlap ---------------------------
c1 <- c()
c2 <- c()
oval <- c()
pval <- c()
pvalu <- c()
loop1 <- length(all_cancer)-1
loop2 <- length(all_cancer)

for(ov1 in 1:loop1){
    m <- ov1+1
    for(ov2 in m:loop2){
        e1 <- length(TF_splicing_events[[ov1]])
        e2 <- length(TF_splicing_events[[ov2]])
        c1 <- c(c1, paste0(all_cancer[ov1],' (', e1,')'))
        c2 <- c(c2, paste0(all_cancer[ov2],' (', e2,')'))
        inter_e <- intersect(TF_splicing_events[[ov1]], TF_splicing_events[[ov2]])
        oval <- c(oval, length(inter_e))

        hyp <- phyper(length(inter_e)-1,e1,length(background_as)-e1,e2,lower.tail = FALSE)
        hyp <- signif(hyp, digits=3)
        pval <- c(pval, hyp)

        ## under representation
        hyp <- phyper(length(inter_e),e1,length(background_as)-e1,e2,lower.tail = TRUE)
        hyp <- signif(hyp, digits=3)
        pvalu <- c(pvalu, hyp)
    # }
    }
}

qval <- p.adjust(pval, 'fdr')
qvalu <- p.adjust(pvalu, 'fdr')
qvalx <- rep('Expected by chance',length(pval))
for(i in 1:length(pval)){
    if(qval[i] < fdr){
        qvalx[i] <- 'Significantly \nmore overlapping'
    }
    if(qvalu[i] < fdr){
        qvalx[i] <- 'Significantly \nless overlapping'
    }
}


pdata <- data.frame(c1=c1, c2=c2, val=oval, pval=qvalx, qval=qval)

cols <- c('#377eb8','#e41a1c')#rev(brewer.pal(3,"Spectral"))
p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pval),colour = "white")+
  theme(legend.text=element_text(size=8))+scale_fill_manual(values=cols,drop=FALSE)
basesize <- 8
p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type", y = "Cancer type") +
  scale_y_discrete() +
  scale_x_discrete()+coord_fixed()+
  guides(fill=guide_legend(title=""))+
  geom_text(data=pdata,aes(y=c2,label=val),size=1.5)+
  theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0.5, colour = "black"))
  # +guides(fill='none')
ggsave(p,filename=paste0(save_dir,'/AS_events_TFs_overlap_genes.png'),width=4, height=2.5, dpi=400)



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

}

pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of TFs with at least one\n perturbed splicing event", limits=c(0,(max(pdata$count))+100)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_overlap_genes.png"),width=7, height=3, dpi=400)



##----- unique splicing events -----------------------------------------------------
num_events_unq <- list()
num_events_unq_counts <- c()

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

    temp <- data.table::fread(all_files[k], sep='\t')
    tempy <- temp[which(temp$as_id %in% temp_unq), ]

}


pdata <- data.frame(cancer=all_cancer, count=num_events_unq_counts)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of TFs with \nperturbed splicing events", limits=c(0,25)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Unique_sig_events_TFs_genes.png"),width=3.5, height=3, dpi=400)



##--- Overlap with oncogenes and cancer suppressor genes ----------------
cancergenes <- data.table::fread('../data/cancerGeneList.tsv')
## keep genes with at least two evidence
wh <- which(cancergenes$`# of occurrence within resources (Column J-P)` > 1)
cancergenes <- cancergenes[wh, ]
oncogenes <- cancergenes[cancergenes$`Is Oncogene` == 'Yes', ]
supgenes <- cancergenes[cancergenes$`Is Tumor Suppressor Gene` == 'Yes', ]
tf_onco <- vector(mode = "list", length = length(all_cancer))
tf_sup <- vector(mode = "list", length = length(all_cancer))

for(k in 1:length(all_cancer)){ ## overlap with the TFs in multiple cancer types

    tf_onco[[k]] <- intersect(oncogenes$`Hugo Symbol`, num_events_combs[[k]])
    tf_sup[[k]] <- intersect(supgenes$`Hugo Symbol`, num_events_combs[[k]])

}



# gene_data <- list(TF_splicing_events, num_events_combs, num_events_unq)
# saveRDS(gene_data, '../data/TF_genes_cancer.rds')



## temp <- readRDS('../data/TF_genes_cancer.rds')

# ###-- which events capture the two TFs affected in all 15 cancer types ---
# tempx <- temp[[2]][[15]]
# all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI.txt', full.names=TRUE))
# for(k in 1:length(tempx)){

#     for(j in 1:length(all_cancer)){
#         tempf <- data.table::fread(all_files[j])
#         wh1 <- which(tempf$FDR < fdr)
#         wh2 <- which(abs(tempf$MEDIAN_DIFF) > diff)
#         wh <- intersect(wh1, wh2)
#         tempfx <- tempf[wh, ]
#         tempfx1 <- tempfx[tempfx$symbol == tempx[k], ]
#     }

# }





