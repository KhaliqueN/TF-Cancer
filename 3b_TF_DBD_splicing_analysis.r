##############################################################################################
# Purpose: For each cancer type, see whether DBDs are removed by AS
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)
# library(dorothea)
library(pheatmap)
library(ggrepel)
library(ggpp)


save_dir <- '../results_new/TF_DBD'
if(dir.exists(save_dir)){
    unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

psi_input <- '../data/PSI_data'
as_input <- '../data/uniprot_Ensembl_Exon_map_DBD_AS'

fdr <- 0.05
atleast_DBD <- 1 ## least number of DBD overlapping amino acids perturbed by a splicing event X to consider the DBD as perturbed by X
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
tf_ensemb_map <- as.data.frame(data.table::fread('../data/TF_ensembl_uniprot.txt', sep='\t'))
tcga_map <- data.table::fread(paste0(psi_input,'/TCGA_SpliceSeq_Gene_Structure.txt'))

all_filesxx <- gtools::mixedsort(list.files(psi_input, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_filesxx <- all_filesxx[-4]
all_cancer <- substr(basename(all_filesxx), 1,4)
paired_sam <- data.table::fread('../data/cancer_paired_samples.txt')
paired_sam <- paired_sam[-4]
## Some of the event coordinates from TCGA splice seq does not overlap to any of the exons mapped from the canonical protein
## sequence from Uniprot. This is because the canonical seqeunce is not always the longest.

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

tevent <- c()
taa <- c()
tcancer <- c()

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
        temp1 <- temp[(temp$ES != '')&(!is.na(temp$ES)), ]
        temp2 <- temp[(temp$AP != '')&(!is.na(temp$AP)), ]
        temp3 <- temp[(temp$AT != '')&(!is.na(temp$AT)), ]
        temp4 <- temp[(temp$AD != '')&(!is.na(temp$AD)), ]
        temp5 <- temp[(temp$AA != '')&(!is.na(temp$AA)), ]
        temp6 <- temp[(temp$ME != '')&(!is.na(temp$ME)), ]

        if(nrow(temp1) > 0){
            wht <- which(temp1$DBD != '-')
            if(length(wht) >= atleast_DBD){
                temp1 <- temp1[wht,]
                all_filesx1 <- c(all_filesx1, all_files[j])
                events1 <- union(events1, unlist(strsplit(as.character(temp1$ES),';')) ) ##--
                ttx <- plyr::count(unlist(strsplit(as.character(temp1$ES),';')) ) ##--
                tevent <- c(tevent, ttx$x)
                taa <- c(taa, ttx$freq)
                tcancer <- c(tcancer, rep(all_cancer[k], length(ttx[[1]])))
                # break
            }
        }

        if(nrow(temp2) > 0){
            wht <- which(temp2$DBD != '-')
            if(length(wht) >= atleast_DBD){
                temp2 <- temp2[wht,]
                all_filesx2 <- c(all_filesx2, all_files[j])
                events2 <- union(events2, unlist(strsplit(as.character(temp2$AP),';')))
                ttx <- plyr::count(unlist(strsplit(as.character(temp2$AP),';')))
                tevent <- c(tevent, ttx$x)
                taa <- c(taa, ttx$freq)
                tcancer <- c(tcancer, rep(all_cancer[k], length(ttx[[1]])))
            }
        }

        if(nrow(temp3) > 0){
            wht <- which(temp3$DBD != '-')
            if(length(wht) >= atleast_DBD){
                temp3 <- temp3[wht,]
                all_filesx3 <- c(all_filesx3, all_files[j])
                events3 <- union(events3, unlist(strsplit(as.character(temp3$AT),';')))
                ttx <- plyr::count(unlist(strsplit(as.character(temp3$AT),';')))
                tevent <- c(tevent, ttx$x)
                taa <- c(taa, ttx$freq)
                tcancer <- c(tcancer, rep(all_cancer[k], length(ttx[[1]])))
            }
        }

        if(nrow(temp4) > 0){
            wht <- which(temp4$DBD != '-')
            if(length(wht) >= atleast_DBD){
                temp4 <- temp4[wht,]
                all_filesx4 <- c(all_filesx4, all_files[j])
                events4 <- union(events4, unlist(strsplit(as.character(temp4$AD),';')))
                ttx <- plyr::count(unlist(strsplit(as.character(temp4$AD),';')))
                tevent <- c(tevent, ttx$x)
                taa <- c(taa, ttx$freq)
                tcancer <- c(tcancer, rep(all_cancer[k], length(ttx[[1]])))
            }
        }

        if(nrow(temp5) > 0){
            wht <- which(temp5$DBD != '-')
            if(length(wht) >= atleast_DBD){
                temp5 <- temp5[wht,]
                all_filesx5 <- c(all_filesx5, all_files[j])
                events5 <- union(events5, unlist(strsplit(as.character(temp5$AA),';')))
                ttx <- plyr::count(unlist(strsplit(as.character(temp5$AA),';')))
                tevent <- c(tevent, ttx$x)
                taa <- c(taa, ttx$freq)
                tcancer <- c(tcancer, rep(all_cancer[k], length(ttx[[1]])))
            }
        }

        if(nrow(temp6) > 0){
            wht <- which(temp6$DBD != '-')
            if(length(wht) >= atleast_DBD){
                temp6 <- temp6[wht,]
                all_filesx6 <- c(all_filesx6, all_files[j])
                events6 <- union(events6, unlist(strsplit(as.character(temp6$ME),';')))
                ttx <- plyr::count(unlist(strsplit(as.character(temp6$ME),';')))
                tevent <- c(tevent, ttx$x)
                taa <- c(taa, ttx$freq)
                tcancer <- c(tcancer, rep(all_cancer[k], length(ttx[[1]])))
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

aa_purt <- data.frame(CANCER=tcancer, ASID=tevent, AA=taa)
aa_purtx <- unique(aa_purt[,c(2,3)])

## check the lengths
# yy = mapply(function(x,y) length(x) > length(y), TFs_with_affected_DBD_ES, ES_affecting_TFs_with_DBD)
##----------------------------------------------------------------------------
##--- which TFs have a DBD info ----
counter <- 0
for(k in 1:length(all_files)){
    temp <- data.table::fread(all_files[k])
    wh <- which(temp$DBD != '-')
    if(length(wh)!=0){
        counter <- counter+1
    }
}
##----------------------------------

# #----- All TFs and all events ---
# all_tfs <- list()
# all_events <- list()
# for(k in 1:length(all_cancer)){
#     temptf <- c()
#     tempev <- c()
#     for(j in 1:length(TFs_list)){
#         temptf <- union(temptf, TFs_list[[j]][[k]])
#         tempev <- union(tempev, Events_list[[j]][[k]])
#     }
#     all_tfs[[k]] <- temptf
#     all_events[[k]] <- tempev
# }

##--------------------------------------------------------------------

##--- Save all DBD perturbing events -----------------------------

##--- Which of the event pairs are also have affected DBDs -----
store_val <- function(y1){
    if(nrow(y1) == 0){
        flag <- 0
    }else{
        y2 <- y1[y1$DBD != '-', ]
        if(nrow(y2) < atleast_DBD){ 
            flag <- 0
        }else{
            flag <- unique(y2$DBD)
        }
    }
    return(paste(flag,collapse=','))
}


odbd <- c()
tcancer <- c()
asid <- c()
corr1 <- c()
corr2 <- c()
corr3 <- c()
corr4 <- c()
counter <- 0
tgenet <- c()
tsplice <- c()
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Events_perturbing_DBDs.xlsx'))

for(k in 1:length(all_cancer)){

    tempf <- data.table::fread(all_filesxx[k], sep='\t')
    ect <- c()
    for(j in 1:length(Events_list)){
        ect <- union(ect, Events_list[[j]][[k]])
    }

    if(length(ect) != 0){

        tempyy <- tempf[tempf$as_id %in% ect, ]

        odbdt <- c()
        asidt <- c()
        corr1t <- c()
        corr2t <- c()
        corr3t <- c()
        corr4t <- c()
        tgene <- c()

        for(j in 1:length(tempyy[[1]])){
        
            tuni <- tf_ensemb_map[tf_ensemb_map$HGNC_symbol == tempyy$symbol[j],]$Uniprotswissprot
            if(length(tuni) == 0){
                counter <- counter+1
                next
            } ## some gene symbols of TFs from TCGA SpliceSeq were not mapped to Ensembl
            tmap <- as.data.frame(data.table::fread(paste0('../data/uniprot_Ensembl_Exon_map_DBD_AS/',tuni,'_',all_cancer[k],'.txt')))
            
            if(tempyy$splice_type[j] == 'AP'){
                tempy1 <- strsplit(as.character(tmap$AP),';')
                why1 <- which(unlist(lapply(tempy1, function(x) tempyy$as_id[j] %in% x)))  
                y1 <- tmap[why1, ]
                odbd <- c(odbd, store_val(y1)) 
                odbdt <- c(odbdt, store_val(y1)) 
            }else if(tempyy$splice_type[j] == 'AT'){
                tempy1 <- strsplit(as.character(tmap$AT),';')
                why1 <- which(unlist(lapply(tempy1, function(x) tempyy$as_id[j] %in% x))) 
                y1 <- tmap[why1, ]
                odbd <- c(odbd, store_val(y1)) 
                odbdt <- c(odbdt, store_val(y1))  
            }else if(tempyy$splice_type[j] == 'ES'){
                tempy1 <- strsplit(as.character(tmap$ES),';')
                why1 <- which(unlist(lapply(tempy1, function(x) tempyy$as_id[j] %in% x))) 
                y1 <- tmap[why1, ]
                odbd <- c(odbd, store_val(y1)) 
                odbdt <- c(odbdt, store_val(y1))   
            }else if(tempyy$splice_type[j] == 'AD'){
                tempy1 <- strsplit(as.character(tmap$AD),';')
                why1 <- which(unlist(lapply(tempy1, function(x) tempyy$as_id[j] %in% x))) 
                y1 <- tmap[why1, ]
                odbd <- c(odbd, store_val(y1)) 
                odbdt <- c(odbdt, store_val(y1))  
            }else if(tempyy$splice_type[j] == 'AA'){
                tempy1 <- strsplit(as.character(tmap$AA),';')
                why1 <- which(unlist(lapply(tempy1, function(x) tempyy$as_id[j] %in% x))) 
                y1 <- tmap[why1, ]
                odbd <- c(odbd, store_val(y1))  
                odbdt <- c(odbdt, store_val(y1)) 
            }else if(tempyy$splice_type[j] == 'ME'){
                tempy1 <- strsplit(as.character(tmap$ME),';')
                why1 <- which(unlist(lapply(tempy1, function(x) tempyy$as_id[j] %in% x))) 
                y1 <- tmap[why1, ]
                odbd <- c(odbd, store_val(y1)) 
                odbdt <- c(odbdt, store_val(y1)) 
            }else{
                odbd <- c(odbd, 0) 
                odbdt <- c(odbdt, 0) 
            }
            
            asid <- c(asid, tempyy$as_id[j])
            tcancer <- c(tcancer, all_cancer[k])
            corr1 <- c(corr1, tempyy$MEAN_CANCER[j])
            corr2 <- c(corr2, tempyy$MEAN_NORMAL[j])
            corr3 <- c(corr3, tempyy$MEAN_DIFF[j])
            corr4 <- c(corr4, tempyy$FDR[j])

            asidt <- c(asidt, tempyy$as_id[j])
            corr1t <- c(corr1t, tempyy$MEAN_CANCER[j])
            corr2t <- c(corr2t, tempyy$MEAN_NORMAL[j])
            corr3t <- c(corr3t, tempyy$MEAN_DIFF[j])
            corr4t <- c(corr4t, tempyy$FDR[j])
            tgene <- c(tgene, tempyy$symbol[j])
            tgenet <- c(tgenet, tempyy$symbol[j])
            tsplice <- c(tsplice, tempyy$splice_type[j])
        }

    }
    
    ## save excel sheet ----
    tdatat <- data.frame(GENE=tgene, AS=asidt, DBD=odbdt, MEAN_CANCER=corr1t, MEAN_NORMAL=corr2t, MEAN_DIFF=corr3t, FDR=corr4t)
    tdatat <- tdatat[tdatat$DBD != 0, ]
    tdatat <- tdatat[order(abs(tdatat$MEAN_DIFF), decreasing=TRUE), ]
    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tdatat)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Events_perturbing_DBDs.xlsx'), overwrite = T)
}

tdata <- data.frame(CANCER=tcancer, SYMBOL=tgenet, AS=asid, SPLICE_TYPE=tsplice, DBD=odbd, MEAN_CANCER=corr1, MEAN_NORMAL=corr2, MEAN_DIFF=corr3, FDR=corr4)
data.table::fwrite(tdata, paste0('../data/Events_perturbing_DBD.txt'), row.names=FALSE, quote=FALSE, sep='\t')

##----------------------------------------------------------------
dbd_purt_dt <- tdata

##--- Events plots ------
sig_events <- list()
tf_events <- list()
as_tags <- c('ES','AP','AT','AD','AA','ME')
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()

for(k in 1:length(all_cancer)){

    tempf <- data.table::fread(all_filesxx[k], sep='\t')
    ect <- c()
    for(j in 1:length(Events_list)){
        ect <- union(ect, dbd_purt_dt[dbd_purt_dt$CANCER == all_cancer[k], ]$AS)
    }

    tempx <- tempf[tempf$as_id %in% ect, ]

    sig_events[[k]] <- tempx$as_id
    tf_events[[k]] <- unique(tempx$symbol)

    temp_count <- plyr::count(tempx$splice_type)
    temp_count$prct <- signif((temp_count$freq/sum(temp_count$freq))*100, 3)
    AA <- c(AA, ifelse(length(which(temp_count$x == 'AA') != 0),temp_count$prct[which(temp_count$x == 'AA')],0))
    AD <- c(AD, ifelse(length(which(temp_count$x == 'AD') != 0),temp_count$prct[which(temp_count$x == 'AD')],0))
    AP <- c(AP, ifelse(length(which(temp_count$x == 'AP') != 0),temp_count$prct[which(temp_count$x == 'AP')],0))
    AT <- c(AT, ifelse(length(which(temp_count$x == 'AT') != 0),temp_count$prct[which(temp_count$x == 'AT')],0))
    ES <- c(ES, ifelse(length(which(temp_count$x == 'ES') != 0),temp_count$prct[which(temp_count$x == 'ES')],0))
    ME <- c(ME, ifelse(length(which(temp_count$x == 'ME') != 0),temp_count$prct[which(temp_count$x == 'ME')],0))
}

##--- plot the number of splicing events affecting TFs -----------------
pdatae <- data.frame(cancer=all_cancer, count1=lengths(sig_events), count2=lengths(tf_events))
maxv <- max(pdatae$count1)
colnames(pdatae) <- c('Cancer','PTSE','TF')
pdata <- reshape2::melt(pdatae)
pdata[pdata == 0] <- NA
p <- ggplot(pdata, aes(Cancer, value, fill=variable)) + 
geom_bar(stat="identity",position="dodge")+
theme(legend.text=element_text(size=12))
basesize <- 8
p <- p + theme_bw(base_size = basesize) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of PTSEs or TFs", limits=c(0,maxv+10)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
scale_fill_manual(values=c('#d95f02','#1b9e77'))+
geom_text(aes(label=value), position=position_dodge(width=0.9),hjust=0, vjust=0.5, angle=85, size=3)+
theme(axis.text.x = element_text(size = basesize, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize), axis.title=element_text(basesize*1.25), legend.position=c(0.85,0.82))+
guides(fill=guide_legend(title="Entity",ncol=1))
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs.png"),width=3.5, height=3, dpi=400)


##--- plot the number of splicing events of different types affecting TFs -----------------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, X=lengths(sig_events), Cancer=all_cancer)
pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","F","X"))  # setting level order
basesize <- 8
ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
# geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
scale_y_continuous(limits=c(0,101), breaks = seq(0, 100, by = 10))+
xlab("Cancer type")+ylab("% of PTSEs affecting DBDs")+
scale_fill_manual(labels=c("A" = "Alternate Acceptor", 
    "B"="Alternate Donor", "C"="Alternate Promoter","D"="Alternate Terminator","E"="Exon skip","F"="Mutually Exclusive Exons", "G"="Retained Intron"), 
values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"), legend.position="bottom",axis.title=element_text(size=basesize),
    legend.text=element_text(size=basesize), legend.title=element_text(size=basesize))+
guides(fill=guide_legend(title="Type of\nalternative\nsplicing event", ncol=1, override.aes = list(size = 1)))
ggsave(ppx,filename=paste0(save_dir, "/Events_types_perturbing_DBDs.png"),width=3, height=4, dpi=400)


###------ Distribution of perturbation values vs. DBD AAs perturbed ----------------------------------------------------
##------------------------------------------------
taa <- c()
tcancer <- c()
tdiff <- c()
tgene <- c()
tsplice <- c()
tdbd <- c()
tas <- c()

for(k in 1:length(all_cancer)){

    temp1 <- dbd_purt_dt[dbd_purt_dt$CANCER == all_cancer[k], ]
    temp2 <- aa_purt[aa_purt$CANCER == all_cancer[k], ]

    for(j in 1:length(temp1[[1]])){
        if(temp1$AS[j] %in% temp2$ASID){
            taa <- c(taa, temp2[temp2$ASID == temp1$AS[j], ]$AA)
            tcancer <- c(tcancer, all_cancer[k])
            tdiff <- c(tdiff, temp1$MEAN_DIFF[j])
            tgene <- c(tgene, temp1$SYMBOL[j])
            tdbd <- c(tdbd, temp1$DBD[j])
            tsplice <- c(tsplice, temp1$SPLICE_TYPE[j])
            tas <- c(tas, temp1$AS[j])
        }
    }
}

pdata <- data.frame(CANCER=tcancer, GENE=tgene, ASID=tas, SPLICE_TYPE=tsplice, DBD=tdbd, AA=taa, MEAN_DIFF=tdiff)
pdata$ID <- paste0(pdata$GENE,'_',pdata$ASID,'_',pdata$SPLICE_TYPE,'\n','(',pdata$CANCER,', ',pdata$DBD,')')

tolabel <- subset(pdata, abs(MEAN_DIFF) > 0.4)$ID
whl <- which(pdata$ID %in% tolabel)
pdata$PL <- ""
pdata$PL[whl] <- tolabel

p <- ggplot(pdata, aes(AA, MEAN_DIFF, color=CANCER, label=PL)) + 
geom_point()+
theme(legend.text=element_text(size=12))
basesize <- 8
p <- p + 
scale_x_continuous(name="# of DBD residues") + 
scale_y_continuous(name="Mean \u0394PSI", limits=c(-0.6,0.6)) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#fb9a99','#ffff99',
    '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
geom_text_repel(position = 
position_nudge_center(x = 0.25,
                      y = 0.25,
                      center_x = 0,
                      center_y = 2,
                      direction = "radial"),
                      family = "Poppins",
                      size = 2.5,
                      color='black',
                      arrow = arrow(length = unit(0.010, "npc")),
                      min.segment.length = 0,
                      hjust = "outward", vjust = "outward") +
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"), legend.position="bottom",axis.title=element_text(size=basesize),
    legend.text=element_text(size=basesize), legend.title=element_text(size=basesize))+
guides(color='none')
ggsave(p,filename=paste0(save_dir,"/AA_vs_meanDiff.png"),width=4, height=3, dpi=500)

# geom_text_repel( family = "Poppins",
#     size = 2,
#     min.segment.length = 0, 
#     seed = 42, 
#     box.padding = 0.5,
#     max.overlaps = Inf,
#     arrow = arrow(length = unit(0.010, "npc")),
#     nudge_x = 0,
#     nudge_y = 0,
#     color = "black") +





# ##--- Splicing events occuring in multiple cancer types --------------------------------
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
#     temp_unn <- c()
#     tempy <- data.frame(matrix(ncol=0, nrow=0))

#     for(i in 1:loop1){
#         temp_ovl <- sig_events[[which(all_cancer == temp_comb[[i]][1])]]
#         if(loop2 > 1){
#             for(j in 2:loop2){
#                 temp_ovl <- intersect(temp_ovl, sig_events[[which(all_cancer == temp_comb[[i]][j])]])
#             }
#         }
#         temp_unn <- union(temp_unn, temp_ovl)
#     }
#     num_events_combs[[k]] <- temp_unn
#     num_events_combs_counts <- c(num_events_combs_counts, length(temp_unn))
# }
# pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)

# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity",position=position_dodge())+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# maxv <- max(pdata$count)
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="# of cancer types") + 
# scale_y_continuous(name="# of splicing events", limits=c(0,maxv+20)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
# ggsave(p,filename=paste0(save_dir,"/Events_perturbing_DBDs_overlap.png"),width=3.5, height=3, dpi=400)



# ##--- plot the number of TFs -----------------
# pdata <- data.frame(cancer=all_cancer, count=lengths(all_tfs))
# p <- ggplot(pdata, aes(cancer, count)) + 
# geom_bar(stat="identity",position=position_dodge())+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="Cancer type") + 
# scale_y_continuous(name="# of TFs with DBDs \naffected by a splicing event", limits=c(0,max(pdata$count)+5)) +
# geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')
# ggsave(p,filename=paste0(save_dir,"/TFs_with_affected_DBDs.png"),width=3.5, height=3, dpi=400)


# ##----- PSI profile of the pancancer (at least 10 cancer types) events -------------
# events_to_consider <- num_events_combs[[10]]
# tcancer <- c()
# tvalue <- c()
# tevent <- c()
# for(k in 1:length(all_cancer)){
#     temp <- data.table::fread(all_filesxx[k], sep='\t')
#     temp1 <- temp[temp$as_id %in% events_to_consider, ]
#     temp1$ID <- paste0(temp1$symbol,'_',temp1$as_id,'_',temp1$splice_type)
#     tcancer <- c(tcancer, rep(all_cancer[k], length(temp1[[1]])))
#     tvalue <- c(tvalue, temp1$MEDIAN_DIFF)
#     tevent <- c(tevent, temp1$ID)
# }
# pdata <- data.frame(tcancer=tcancer, val=tvalue, event=tevent)

# pdata_mat <- as.data.frame(matrix(nrow=length(unique(pdata$event)), ncol=length(all_cancer),0))
# rownames(pdata_mat) <- gtools::mixedsort(unique(pdata$event))
# colnames(pdata_mat) <- gtools::mixedsort(unique(all_cancer))

# for(k in 1:length(pdata[[1]])){
#     wh1 <- which(rownames(pdata_mat) == pdata$event[k])
#     wh2 <- which(colnames(pdata_mat) == pdata$tcancer[k])
#     pdata_mat[wh1, wh2] <- as.numeric(pdata$val[k])
# }

# # pdata_mat <- cbind(data.frame(Event=unique(pdata$event)), pdata_mat)
# pdx <- t(as.matrix(pdata_mat))

# color <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
# breaks <- seq(-1, 1, 0.2)
# p <- pheatmap(pdx,fontsize=3, color=color, breaks=breaks, cluster_rows=FALSE, cluster_cols=FALSE,cellheight=5, cellwidth = 5)
# ggsave(p,filename=paste0(save_dir, "/Pancancer_events_DBD.png"),width=7, height=5, dpi=600)




##--- total number of events -----
# tevs <- unique(unlist(sig_events)) ## 316

##--- Distribution of DBD types -----------------------------------------
all_files <- list.files(as_input, pattern=paste0('*',all_cancer[1],'.txt'), full.names=TRUE)
all_dbds <- c()
numdbd <- c()
# rel_pos <- list()
# counter <- 1
for(k in 1:length(all_files)){
    temp <- data.table::fread(all_files[k])
    tempdb <- setdiff(unique(temp$DBD),'-')
    if(length(tempdb) != 0){
        wh <- which(tempdb == 'NR C4-type')
        if(length(wh) != 0){
            tempdb[wh] <- 'Nuclear receptor'
        }
    }
    tempdb <- unique(tempdb)
    all_dbds <- c(all_dbds, tempdb)
    numdbd <- c(numdbd, length(tempdb))
    # wh <- which(temp$DBD != '-')
    # if(length(wh) != 0){
    #     rel_pos[[counter]] <- ceiling(signif(wh/length(temp[[1]]), 1))
    #     counter <- counter+1
    # }
}

# ## which TFs have more than one DBD --> check them manually
# whm <- which(numdbd > 1)
# for(k in 11:15){
#     temp <- data.table::fread(all_files[whm[k]])
#     print(temp$DBD)
# }

pdata <- plyr::count(all_dbds)
pdata <- pdata[order(-pdata$freq),]
wh <- which(pdata$freq < 5)
pdata2 <- pdata[wh,]
pdata1 <- pdata[-wh,]
pdatax <- rbind(pdata1, data.frame(x=paste0('Others (',length(pdata2[[1]]),')'),freq=sum(pdata2$freq)))
pdatax$x <- factor(pdatax$x, levels=pdatax$x)

p <- ggplot(pdatax, aes(x, freq)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="DNA binding domain") + 
scale_y_continuous(name="Frequency of occurence", limits=c(0,max(pdatax$freq)+40)) +
geom_text(aes(label=freq), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 1,vjust=1, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/DBDs_background.png"),width=7, height=4, dpi=400)





##--- Distribution of DBD positions --------------------------------------------------------------------
rel_pos <- c()
counter <- 1
idbd <- c()
ids <- c()
for(k in 1:length(all_files)){
    temp <- data.table::fread(all_files[k])
    tempdb <- setdiff(unique(temp$DBD),'-')
    if(length(tempdb) != 0){
        wh <- which(tempdb == 'NR C4-type')
        if(length(wh) != 0){
            tempdb[wh] <- 'Nuclear receptor'
        }
    }
    tempdb <- unique(tempdb)

    wh <- which(temp$DBD != '-')

    if(length(wh) != 0){
        whi <- intersect(tempdb, all_dbds)
        for(j in 1:length(whi)){
            whx <- which(temp$DBD == whi[j])
            rel_pos <- c(rel_pos,whx/length(temp[[1]]))
            counter <- counter+1
            idbd <- c(idbd, rep(whi[j], length(whx)))
            ids <- c(ids, rep(counter, length(whx)))
        }
    }
}

qdata <- data.frame(DBD=idbd, RPOS=rel_pos, ID=ids)
## bin the relative positions
qdata$bval <- cut(qdata$RPOS, breaks = seq(0,1,0.05), include.lowest=TRUE)

wh <- which(qdata$DBD %in% pdata1$x)
whe <- setdiff(seq(1,length(qdata[[1]])), wh)
qdata$DBD[whe] <- paste0('Others (', length(pdata2[[1]]), ')')

udbd <- unique(qdata$bval)
qqdata <- data.frame(matrix(nrow=0, ncol=4))
for(j in 1:length(udbd)){
    temp <- qdata[qdata$bval == udbd[j], ]
    tempu <- plyr::count(temp$DBD)
    tempu$PRT <- (tempu$freq/sum(tempu$freq))*100
    tempu$FLAG <- rep(udbd[j],length(tempu[[1]]))
    qqdata <- rbind(qqdata, tempu)
}

basesize <- 8
qqdata$lfreq <- log2(qqdata$freq)
qqdata$x <- factor(qqdata$x, levels=pdatax$x)

cols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
    '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#000000',
    '#ffffff','#636363')
p <- ggplot(qqdata, aes(FLAG, freq, fill=x)) + geom_bar(stat="identity",position="stack",color='black')
p <- p + theme_bw() +
  scale_y_continuous(name="# of amino acids") +
  scale_x_discrete(name="Relative sequence position bracket") +
  scale_fill_manual(values=cols)+
  guides(fill=guide_legend(title="DNA binding domain type",ncol=2))+
  theme(axis.text.x = element_text(size = basesize * 1,angle = 60, hjust = 1,vjust=1, colour = "black"), 
                         axis.text.y = element_text(size = basesize * 1, colour='black'), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=basesize*1),axis.title.y = element_text(size=basesize*1), 
        plot.title = element_text(size=basesize*1), strip.text.x = element_text(size = basesize * 1, colour = "black", angle = 0), 
        strip.text.y = element_text(size = basesize * 1, colour = "black", angle = 0), legend.position="bottom", 
        legend.text=element_text(size=10))+guides(fill='none')
ggsave(p,filename=paste0(save_dir,'/DBD_distribution.png'),width=3.5, height=4, dpi=500)

# ## distribution of amino acid positions ---
# basesize <- 8
# p <- ggplot(qdata, aes(RPOS)) + geom_barplot(binwidth=0.01)
# p <- p + theme_grey(base_size = basesize) + labs(x = "Relative protein sequence position", y = "DBD type") +
#   scale_y_continuous(name="DBD type") +
#   scale_x_continuous(name="Relative protein sequence position") +
#   guides(fill='none')+
#   theme(axis.text.x = element_text(size = basesize * 1,angle = 60, hjust = 1,vjust=1, colour = "black"), 
#                          axis.text.y = element_text(size = basesize * 1, colour='black'), 
#                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.title.x = element_text(size=basesize*1),axis.title.y = element_text(size=basesize*1), 
#         plot.title = element_text(size=basesize*1), strip.text.x = element_text(size = basesize * 1, colour = "black", angle = 0), 
#         strip.text.y = element_text(size = basesize * 1, colour = "black", angle = 0), legend.text=element_text(size=10))
# ggsave(p,filename=paste0(save_dir,'/DBD_distribution.png'),width=5, height=5, dpi=500)



##-----------------------------------

##--- which DBD types are affected over all cancer --------------------------------------
##---------------------------------------------------------------------------------------
##-- unique event DBD pair --
event_dbd_unq <- unique(dbd_purt_dt[, c(3,5)])
wh <- which(event_dbd_unq$DBD %like% 'NR C4-type')
event_dbd_unq$DBD[wh] <- 'Nuclear receptor'
wh <- which(event_dbd_unq$DBD %like% ',')
pdatamul <- event_dbd_unq[wh,]
pdatauni <- event_dbd_unq[-wh,]
for(k in 1:length(pdatamul[[1]])){
    tempv <- unlist(strsplit(pdatamul$DBD[k], '[,]'))
    for(j in 1:length(tempv)){
        pdatauni <- rbind(pdatauni, data.frame(AS=pdatamul$AS[k], DBD=tempv[j]))
    }
}

pdataunix <- plyr::count(pdatauni$DBD)

wh <- which(pdataunix$x %in% pdata1$x)
whe <- setdiff(seq(1,length(pdataunix[[1]])), wh)
pdataunix$x[whe] <- paste0('Others (', length(pdata2[[1]]), ')')
pdatauniy <- aggregate(pdataunix$freq, by=list(Category=pdataunix$x), FUN=sum)

## add missing DBDs ---
mdbd <- setdiff(pdatax$x, pdatauniy$Category)
pdatauniy <- rbind(pdatauniy, data.frame(Category=mdbd, x=0))
pdatauniy$Category <- factor(pdatauniy$Category, levels=pdatax$x)

p <- ggplot(pdatauniy, aes(Category, x)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="DNA binding domain") + 
scale_y_continuous(name="# of splicing events", limits=c(0,max(pdatauniy$x)+40)) +
geom_text(aes(label=x), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 1,vjust=1, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/DBDs_affected.png"),width=7, height=4, dpi=400)



pdatagh <- data.frame(matrix(ncol=3, nrow=0))
for(k in 1:length(all_cancer)){

    tempdbd <- unique(dbd_purt_dt[dbd_purt_dt$CANCER == all_cancer[k],][,c(3,5)])
    wh <- which(tempdbd$DBD %like% 'NR C4-type')
    if(length(wh) != 0){tempdbd$DBD[wh] <- 'Nuclear receptor'}
    wh <- which(tempdbd$DBD %like% ',')
    if(length(wh) != 0){
        pdatam <- tempdbd[wh,]
        pdatau <- tempdbd[-wh,]
        for(i in 1:length(pdatam[[1]])){
            tempv <- unlist(strsplit(pdatam$DBD[i], '[,]'))
            for(j in 1:length(tempv)){
                pdatau <- rbind(pdatau, data.frame(AS=pdatam$AS[i], DBD=tempv[j]))
            }
        }
    }else{
        pdatau <- tempdbd
    }
    pdataux <- plyr::count(pdatau$DBD)
    pdataux$CANCER <- rep(all_cancer[k], length(pdataux[[1]]))
    pdatagh <- rbind(pdatagh, pdataux)
}


##make the group "others (29)"
wh <- which(pdatagh$x %in% pdata1$x)
whe <- setdiff(seq(1,length(pdatagh[[1]])), wh)
pdatagh$x[whe] <- paste0('Others (', length(pdata2[[1]]), ')')


##-- add number of unique splicing events ---
nnam <- c()
for(k in 1:length(pdatagh[[1]])){
    wh <- which(pdatauniy$Category == pdatagh$x[k])
    nnam <- c(nnam, paste0(pdatauniy$Category[wh],':',pdatauniy$x[wh]))
}
pdatagh$x <- nnam

nnam <- c()
for(k in 1:length(pdatagh[[1]])){
    wh <- which(pdatae$Cancer == pdatagh$CANCER[k])
    nnam <- c(nnam, paste0(pdatae$Cancer[wh],':',pdatae$PTSE[wh]))
}
pdatagh$CANCER <- nnam

all_cancerx <- gtools::mixedsort(unique(pdatagh$CANCER))


##-- Also make scatter plot of background and forground affected DBD --
##--- tile plot -------------------------------------------------------
pdataghx <- aggregate(pdatagh$freq, by=list(pdatagh$x, pdatagh$CANCER), FUN=sum)
colnames(pdataghx) <- c('x','CANCER','freq')
alldbds <- unique(pdataghx$x)
tdbd <- c()
tcancer <- c()
tfrac <- c()
for(k in 1:length(all_cancerx)){
    temp <- pdataghx[pdataghx$CANCER == all_cancerx[k], ]
    temp$per <- temp$freq/sum(temp$freq)
    for(j in 1:length(alldbds)){
        tempx <- temp[temp$x == alldbds[j], ]
        if(nrow(tempx) != 0){
            tfrac <- c(tfrac, tempx$freq)
        }else{
            tfrac <- c(tfrac, 0)
        }
        tcancer <- c(tcancer, all_cancerx[k])
        tdbd <- c(tdbd, alldbds[j])
    }
}
pdatat <- data.frame(CANCER=tcancer, DBD=tdbd, FRAC=tfrac)
pdatat$bval <- ifelse(pdatat$FRAC < 10, pdatat$FRAC, '>=10')
# pdatat$DBD <- factor(pdatat$DBD, levels=pdatax$x)


cols <- c('#800026','#000000','#ffffff','#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026')
# c('#08306b', '#6baed6','#deebf7', '#ffffcc','#fed976','#fd8d3c','#e31a1c','#800026')#'#2171b5', 
# cols <- rev(brewer.pal(11,"Spectral"))
basesize <- 8
p <- ggplot(pdatat, aes(DBD, CANCER)) + geom_tile(aes(fill = bval), color='black')+scale_fill_manual(values=cols)
  # scale_fill_gradientn(colors=cols)
p <- p + theme_grey(base_size = basesize) + labs(x = "Sample", y = "Gene") +
  scale_y_discrete(name="Cancer type") +
  scale_x_discrete(name="DNA binding domain") +
  guides(fill=guide_legend(title="# of PTSEs", size=10, ncol=2, override.aes = list(size = 2)))+
  theme(axis.text.x = element_text(size = basesize * 1,angle = 60, hjust = 1,vjust=1, colour = "black"), 
                         axis.text.y = element_text(size = basesize * 1, colour='black'), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=basesize*1),axis.title.y = element_text(size=basesize*1), 
        plot.title = element_text(size=basesize*1), strip.text.x = element_text(size = basesize * 1, colour = "black", angle = 0), 
        strip.text.y = element_text(size = basesize * 1, colour = "black", angle = 0), 
        legend.text=element_text(size=10))
ggsave(p,filename=paste0(save_dir,'/Perturbed_DBDs_tile.png'),width=5, height=4, dpi=500)



##---- alternative visualization of the fraction of perturbed DBD types ----
pdatax1 <- pdatax
pdatauniy1 <- pdatauniy
pdatax1$frac <- signif((pdatax1$freq/sum(pdatax1$freq)*100),3)
pdatauniy1$frac <- signif((pdatauniy1$x/sum(pdatauniy1$x)*100),3)
colnames(pdatauniy1) <- c('x','freq','frac')
pdatauniy1$FLAG <- rep('Cancer',length(pdatauniy1[[1]]))
pdatax1$FLAG <- rep('Background',length(pdatax1[[1]]))
alldata <- rbind(pdatax1, pdatauniy1)

cols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
    '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#000000',
    '#ffffff','#636363')
p <- ggplot(alldata, aes(FLAG, frac, fill=x)) + 
geom_bar(stat="identity",position='stack', color='black')+
theme(legend.text=element_text(size=12))
basesize <- 8
p <- p + theme_bw() +
scale_x_discrete(name="") + 
scale_y_continuous(name="% of occurence/perturbation") +
scale_fill_manual(values=cols)+
# geom_text(aes(label=x), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
  theme(axis.text.x = element_text(size = basesize * 1,angle = 60, hjust = 1,vjust=1, colour = "black"), 
                         axis.text.y = element_text(size = basesize * 1, colour='black'), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=basesize*1),axis.title.y = element_text(size=basesize*1), 
        plot.title = element_text(size=basesize*1), strip.text.x = element_text(size = basesize * 1, colour = "black", angle = 0), 
        strip.text.y = element_text(size = basesize * 1, colour = "black", angle = 0), 
        legend.text=element_text(size=basesize), legend.title=element_text(size=basesize*1.25))+
guides(fill=guide_legend(title="DNA binding domain type"),ncol=2)
ggsave(p,filename=paste0(save_dir,"/DBDs_affected_frac.png"),width=5, height=5, dpi=500)


##--------------------------------------------------------------------------