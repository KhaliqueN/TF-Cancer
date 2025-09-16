##############################################################################################
# Purpose: For each cancer type, see whether EDs are removed by AS
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

save_dir <- '../results_new/TF_ED'
if(dir.exists(save_dir)){
    unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

psi_input <- '../data/PSI_data'
as_input <- '../data/uniprot_Ensembl_Exon_map_DBD_ED_AS'

atleast_DBD <- 1 ## least number of ED overlapping amino acids perturbed by a splicing event X to consider the DBD as perturbed by X
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

tevents <- c()
# tf_fam <- c()
taa <- c()
tdbd <- c()
tcancer <- c()
tgene <- c()

tempfun <- function(tempo_data, tempxz, xxp){

    # print(as.numeric(xxp))
    temp1_c <- plyr::count(tempo_data[,c(21,as.numeric(xxp))])
    colnames(temp1_c) <- c('ED','AS','freq')

    xdbd <- c()
    for(ii in 1:length(temp1_c[[1]])){
        xdbd <- c(xdbd, length(which(tempxz$ED == temp1_c$ED[ii])))
    } 
    temp1_c$len <- xdbd
    temp1_c$frac <- signif((temp1_c$freq/temp1_c$len)*100,3)

    wh <- which(temp1_c[[2]] %like% ';')
    if(length(wh) != 0){
        temp1_c1 <- temp1_c[wh,]
        temp1_c2 <- temp1_c[-wh,][,c(1,2,5)]
        for(i in 1:length(temp1_c1[[1]])){
            tempv <- unlist(strsplit(temp1_c1[[2]][i], '[;]'))
            for(j in 1:length(tempv)){
                temp1_c2 <- rbind(temp1_c2, data.frame(ED=temp1_c1$ED[i], AS=tempv[j], frac=temp1_c1$frac[i]))
            }
        }
    }else{
        temp1_c2 <- temp1_c[,c(1,2,5)]
    }

    temp1_cc2 <- aggregate(temp1_c2$frac, by=list(ED=temp1_c2$ED, AS=temp1_c2$AS), FUN=sum)
    colnames(temp1_cc2) <- c('ED','AS','frac')
    return(temp1_cc2)

}

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
            wht <- which(temp1$ED != '-')
            if(length(wht) >= atleast_DBD){
                temp1 <- temp1[wht,]
                all_filesx1 <- c(all_filesx1, all_files[j])
                events1 <- union(events1, unlist(strsplit(as.character(temp1$ES),';')) ) ##--
                tempxg <- tempfun(as.data.frame(temp1), as.data.frame(temp), 15)
                tdbd <- c(tdbd, as.character(tempxg$ED))
                tevents <- c(tevents, as.character(tempxg$AS))
                taa <- c(taa, as.character(tempxg$frac))
                tcancer <- c(tcancer, rep(all_cancer[k], length(tempxg[[1]])))
                tgene <- c(tgene, rep(strsplit(basename(all_files[j]),'[_]')[[1]][1], length(tempxg[[1]])))
            }
        }

        if(nrow(temp2) > 0){
            wht <- which(temp2$ED != '-')
            if(length(wht) >= atleast_DBD){
                temp2 <- temp2[wht,]
                all_filesx2 <- c(all_filesx2, all_files[j])
                events2 <- union(events2, unlist(strsplit(as.character(temp2$AP),';')))
                tempxg <- tempfun(as.data.frame(temp2), as.data.frame(temp), 16)
                tdbd <- c(tdbd, as.character(tempxg$ED))
                tevents <- c(tevents, as.character(tempxg$AS))
                taa <- c(taa, as.character(tempxg$frac))
                tcancer <- c(tcancer, rep(all_cancer[k], length(tempxg[[1]])))
                tgene <- c(tgene, rep(strsplit(basename(all_files[j]),'[_]')[[1]][1], length(tempxg[[1]])))
            }
        }

        if(nrow(temp3) > 0){
            wht <- which(temp3$ED != '-')
            if(length(wht) >= atleast_DBD){
                temp3 <- temp3[wht,]
                all_filesx3 <- c(all_filesx3, all_files[j])
                events3 <- union(events3, unlist(strsplit(as.character(temp3$AT),';')))
                tempxg <- tempfun(as.data.frame(temp3), as.data.frame(temp), 17)
                tdbd <- c(tdbd, as.character(tempxg$ED))
                tevents <- c(tevents, as.character(tempxg$AS))
                taa <- c(taa, as.character(tempxg$frac))
                tcancer <- c(tcancer, rep(all_cancer[k], length(tempxg[[1]])))
                tgene <- c(tgene, rep(strsplit(basename(all_files[j]),'[_]')[[1]][1], length(tempxg[[1]])))
            }
        }

        if(nrow(temp4) > 0){
            wht <- which(temp4$ED != '-')
            if(length(wht) >= atleast_DBD){
                temp4 <- temp4[wht,]
                all_filesx4 <- c(all_filesx4, all_files[j])
                events4 <- union(events4, unlist(strsplit(as.character(temp4$AD),';')))
                tempxg <- tempfun(as.data.frame(temp4), as.data.frame(temp), 18)
                tdbd <- c(tdbd, as.character(tempxg$ED))
                tevents <- c(tevents, as.character(tempxg$AS))
                taa <- c(taa, as.character(tempxg$frac))
                tcancer <- c(tcancer, rep(all_cancer[k], length(tempxg[[1]])))
                tgene <- c(tgene, rep(strsplit(basename(all_files[j]),'[_]')[[1]][1], length(tempxg[[1]])))
            }
        }

        if(nrow(temp5) > 0){
            wht <- which(temp5$ED != '-')
            if(length(wht) >= atleast_DBD){
                temp5 <- temp5[wht,]
                all_filesx5 <- c(all_filesx5, all_files[j])
                events5 <- union(events5, unlist(strsplit(as.character(temp5$AA),';')))
                tempxg <- tempfun(as.data.frame(temp5), as.data.frame(temp), 19)
                tdbd <- c(tdbd, as.character(tempxg$ED))
                tevents <- c(tevents, as.character(tempxg$AS))
                taa <- c(taa, as.character(tempxg$frac))
                tcancer <- c(tcancer, rep(all_cancer[k], length(tempxg[[1]])))
                tgene <- c(tgene, rep(strsplit(basename(all_files[j]),'[_]')[[1]][1], length(tempxg[[1]])))
            }
        }

        if(nrow(temp6) > 0){
            wht <- which(temp6$ED != '-')
            if(length(wht) >= atleast_DBD){
                temp6 <- temp6[wht,]
                all_filesx6 <- c(all_filesx6, all_files[j])
                events6 <- union(events6, unlist(strsplit(as.character(temp6$ME),';')))
                tempxg <- tempfun(as.data.frame(temp6), as.data.frame(temp), 20)
                tdbd <- c(tdbd, as.character(tempxg$ED))
                tevents <- c(tevents, as.character(tempxg$AS))
                taa <- c(taa, as.character(tempxg$frac))
                tcancer <- c(tcancer, rep(all_cancer[k], length(tempxg[[1]])))
                tgene <- c(tgene, rep(strsplit(basename(all_files[j]),'[_]')[[1]][1], length(tempxg[[1]])))
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

aa_purt <- data.frame(CANCER=tcancer, ASID=tevents, AA=taa, DBD=tdbd, GENE=tgene)



# ##--- transfer DBD names----
# temp_dbd_map <- data.table::fread('../data/Prosite_DBD_map.txt')
# mdbd <- c()
# for(k in 1:length(aa_purt[[1]])){
#     tempdd <- paste(setdiff(unique(temp_dbd_map[temp_dbd_map$DBD == aa_purt$DBD[k], ]$PROSITE),''),collapse=';')
#     if(length(tempdd) > 1){break}
#     mdbd <- c(mdbd, tempdd)
# }

# ##--------------------------
# aa_purt$DBD <- unlist(substr(aa_purt$DBD,1,8))
# aa_purtx <- unique(aa_purt[,c(2,3)])

## check the lengths
# yy = mapply(function(x,y) length(x) > length(y), TFs_with_affected_DBD_ES, ES_affecting_TFs_with_DBD)
##----------------------------------------------------------------------------
##--- which TFs have a DBD info ----
counter <- 0
for(k in 1:length(all_files)){
    temp <- data.table::fread(all_files[k])
    wh <- which(temp$ED != '-')
    if(length(wh)!=0){
        counter <- counter+1
    }
}
##----------------------------------

##--------------------------------------------------------------------

##--- Save all DBD perturbing events -----------------------------

##--- Which of the event pairs are also have affected DBDs -----
store_val <- function(y1){
    if(nrow(y1) == 0){
        flag <- 0
    }else{
        y2 <- y1[y1$ED != '-', ]
        if(nrow(y2) < atleast_DBD){ 
            flag <- 0
        }else{
            flag <- unique(y2$ED)
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
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/PTSEs_perturbing_EDs.xlsx'))

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
            tmap <- as.data.frame(data.table::fread(paste0('../data/uniprot_Ensembl_Exon_map_DBD_ED_AS/',tuni,'_',all_cancer[k],'.txt')))
            
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
    # tdatat <- data.frame(GENE=tgene, AS=asidt, DBD=odbdt, MEAN_CANCER=corr1t, MEAN_NORMAL=corr2t, MEAN_DIFF=corr3t, FDR=corr4t)
    tdatat <- data.frame(AS_ID=asidt, DBD=odbdt)
    tdatat <- tdatat[tdatat$DBD != 0, ]
    # tdatat <- tdatat[order(abs(tdatat$MEAN_DIFF), decreasing=TRUE), ]
    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tdatat)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/PTSEs_perturbing_EDs.xlsx'), overwrite = T)
}

tdata <- data.frame(CANCER=tcancer, SYMBOL=tgenet, AS=asid, SPLICE_TYPE=tsplice, ED=odbd, MEAN_CANCER=corr1, MEAN_NORMAL=corr2, MEAN_DIFF=corr3, FDR=corr4)
data.table::fwrite(tdata, paste0('../data/Events_perturbing_ED.txt'), row.names=FALSE, quote=FALSE, sep='\t')
## tempd <- data.table::fread('../data/Events_perturbing_ED.txt')

##----------------------------------------------------------------
dbd_purt_dt <- tdata

##--overlap of AS 
# tempovr <- setdiff(aa_purt$ASID, tdata$AS) ## only one event is missed due to missed mapping of uniprot and ensembl

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
strip.text = element_text(size = basesize), axis.title=element_text(size=basesize*1.25), legend.position=c(0.85,0.82))+
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
ggsave(ppx,filename=paste0(save_dir, "/Events_types_perturbing_EDs.png"),width=3, height=4, dpi=400)


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

    for(j in 1:length(temp2[[1]])){
        if(temp2$ASID[j] %in% temp1$AS){
            taa <- c(taa, temp2$AA[j])
            tcancer <- c(tcancer, temp2$CANCER[j])
            tdiff <- c(tdiff, temp1[temp1$AS == temp2$ASID[j], ]$MEAN_DIFF)
            tgene <- c(tgene, temp1[temp1$AS == temp2$ASID[j], ]$SYMBOL)
            tdbd <- c(tdbd, temp2$DBD[j])
            tsplice <- c(tsplice, temp1[temp1$AS == temp2$ASID[j], ]$SPLICE_TYPE)
            tas <- c(tas, temp2$ASID[j])
        }
    }
}

pdata <- data.frame(CANCER=tcancer, GENE=tgene, ASID=tas, SPLICE_TYPE=tsplice, DBD=tdbd, AA=taa, MEAN_DIFF=tdiff)

pdata$ID <- paste0(pdata$GENE,'_',pdata$ASID,'_',pdata$SPLICE_TYPE,'\n','(',pdata$CANCER,', ',pdata$DBD,')')
pdata$AA <- as.numeric(pdata$AA)
tolabel <- subset(pdata, abs(MEAN_DIFF) > 0.4)$ID
# tolabel2 <- subset(pdata, AA > 80)$ID
# tolabel <- intersect(tolabel1, tolabel2)
whl <- which(pdata$ID %in% tolabel)
pdata$PL <- ""
pdata$PL[whl] <- tolabel

pdata$AA <- as.numeric(pdata$AA)
# p <- ggplot(pdata, aes(AA, MEAN_DIFF, color=CANCER, label=PL)) + 
p <- ggplot(pdata, aes(AA, MEAN_DIFF, color=CANCER)) + 
geom_point()+
theme(legend.text=element_text(size=12))
basesize <- 10
p <- p + 
scale_x_continuous(name="% of ED residues", limits=c(0,110), breaks = seq(0, 100, by = 20)) + 
scale_y_continuous(name="Mean \u0394PSI", limits=c(-0.8, 0.8)) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#fb9a99','#ffff99',
    '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
# geom_text_repel(family = "Poppins",
#     max.overlaps=Inf,
#                       size = 2,
#                       color='black',
#                       arrow = arrow(length = unit(0.010, "npc")),
#                       min.segment.length = 0) +
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),axis.title=element_text(size=basesize),
    legend.text=element_text(size=basesize), legend.title=element_text(size=basesize),
    panel.border = element_blank())+
# guides(color='none')
guides(color=guide_legend(title="Cancer type", ncol=2))
ggsave(p,filename=paste0(save_dir,"/AA_vs_meanDiff.png"),width=5, height=3, dpi=500)


##--- total number of events -----
# tevs <- unique(unlist(sig_events)) 

##--- Distribution of DBD types -----------------------------------------
all_files <- list.files(as_input, pattern=paste0('*',all_cancer[1],'.txt'), full.names=TRUE)
all_dbds <- c()

for(k in 1:length(all_files)){
    temp <- data.table::fread(all_files[k])
    temp1 <- temp[temp$ED != '-', ]
    if(nrow(temp1) != 0){
        tempdb <- unique(temp1$ED)
        all_dbds <- c(all_dbds, tempdb)
    }
}


pdata <- plyr::count(all_dbds)


# pdata <- pdata[order(-pdata$freq),]
# wh <- which(pdata$freq < 5)
# pdata2 <- pdata[wh,]
# pdata1 <- pdata[-wh,]
# pdatax <- rbind(pdata1, data.frame(x=paste0('Others (',length(pdata2[[1]]),')'),freq=sum(pdata2$freq)))
# pdatax$x <- factor(pdatax$x, levels=pdatax$x)

##--affected DBD per gene ---
event_dbd_unq <- unique(dbd_purt_dt[, c(2,5)])
wh <- which(event_dbd_unq$ED %like% ',')
pdatamul <- event_dbd_unq[wh,]
pdatauni <- event_dbd_unq[-wh,]
for(k in 1:length(pdatamul[[1]])){
    tempv <- unlist(strsplit(pdatamul$ED[k], '[,]'))
    for(j in 1:length(tempv)){
        pdatauni <- rbind(pdatauni, data.frame(SYMBOL=pdatamul$SYMBOL[k], ED=tempv[j]))
    }
}
pdataq <- plyr::count(pdatauni[[2]])
# pdataq <- rbind(pdataq, data.frame(x=setdiff(pdata$x, pdataq$x), freq=0))

pdata <- pdata[order(pdata$x),]
pdataq <- pdataq[order(pdataq$x),]
pdataq$frac <- pdataq$freq/pdata$freq
pdataq <- pdataq[order(pdataq$frac), ]

# p <- ggplot(pdatax, aes(x, freq)) + 
# geom_bar(stat="identity",position=position_dodge())+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="DNA binding domain") + 
# scale_y_continuous(name="Frequency of occurence", limits=c(0,max(pdatax$freq)+40)) +
# geom_text(aes(label=freq), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 1,vjust=1, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')
# ggsave(p,filename=paste0(save_dir,"/DBDs_background.png"),width=7, height=4, dpi=400)





##--- Distribution of ED positions --------------------------------------------------------------------
rel_pos <- c()
counter <- 1
idbd <- c()
ids <- c()
for(k in 1:length(all_files)){
    temp <- data.table::fread(all_files[k])
    tempdb <- setdiff(unique(temp$ED),'-')

    tempdb <- unique(tempdb)

    wh <- which(temp$ED != '-')

    if(length(wh) != 0){
        whi <- intersect(tempdb, all_dbds)
        for(j in 1:length(whi)){
            whx <- which(temp$ED == whi[j])
            rel_pos <- c(rel_pos,whx/length(temp[[1]]))
            counter <- counter+1
            idbd <- c(idbd, rep(whi[j], length(whx)))
            ids <- c(ids, rep(counter, length(whx)))
        }
    }
}

qdata <- data.frame(ED=idbd, RPOS=rel_pos, ID=ids)

## bin the relative positions
qdata$bval <- cut(qdata$RPOS, breaks = seq(0,1,0.05), include.lowest=TRUE)

# wh <- which(qdata$ED %in% pdata$x)
# whe <- setdiff(seq(1,length(qdata[[1]])), wh)
# qdata$DBD[whe] <- paste0('Others (', length(pdata2[[1]]), ')')

udbd <- unique(qdata$bval)
qqdata <- data.frame(matrix(nrow=0, ncol=4))
for(j in 1:length(udbd)){
    temp <- qdata[qdata$bval == udbd[j], ]
    tempu <- plyr::count(temp$ED)
    tempu$PRT <- (tempu$freq/sum(tempu$freq))*100
    tempu$FLAG <- rep(udbd[j],length(tempu[[1]]))
    qqdata <- rbind(qqdata, tempu)
}

basesize <- 8
qqdata$lfreq <- log2(qqdata$freq)
# qqdata$x <- factor(qqdata$x, levels=pdatax$x)

cols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
    '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#000000',
    '#ffffff','#636363','#bdbdbd')
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
ggsave(p,filename=paste0(save_dir,'/ED_distribution.png'),width=3.5, height=4, dpi=500)


# ##-----------------------------------

##--- which DBD types are affected over all cancer --------------------------------------
##---------------------------------------------------------------------------------------
##-- unique event DBD pair --
event_dbd_unq <- unique(dbd_purt_dt[, c(3,5)])
wh <- which(event_dbd_unq$ED %like% ',')
pdatamul <- event_dbd_unq[wh,]
pdatauni <- event_dbd_unq[-wh,]
for(k in 1:length(pdatamul[[1]])){
    tempv <- unlist(strsplit(pdatamul$ED[k], '[,]'))
    for(j in 1:length(tempv)){
        pdatauni <- rbind(pdatauni, data.frame(AS=pdatamul$AS[k], ED=tempv[j]))
    }
}

pdataunix <- plyr::count(pdatauni$ED)
pdatauniy <- pdataunix


# wh <- which(pdataunix$x %in% pdata1$x)
# whe <- setdiff(seq(1,length(pdataunix[[1]])), wh)
# pdataunix$x[whe] <- paste0('Others (', length(pdata2[[1]]), ')')
# pdatauniy <- aggregate(pdataunix$freq, by=list(Category=pdataunix$x), FUN=sum)

# ## add missing DBDs ---
# mdbd <- setdiff(pdatax$x, pdatauniy$Category)
# pdatauniy <- rbind(pdatauniy, data.frame(Category=mdbd, x=0))
# pdatauniy$Category <- factor(pdatauniy$Category, levels=pdatax$x)

# p <- ggplot(pdatauniy, aes(Category, x)) + 
# geom_bar(stat="identity",position=position_dodge())+
# theme(legend.text=element_text(size=12))
# basesize <- 12
# p <- p + theme_bw(base_size = basesize * 0.8) +
# scale_x_discrete(name="DNA binding domain") + 
# scale_y_continuous(name="# of splicing events", limits=c(0,max(pdatauniy$x)+40)) +
# geom_text(aes(label=x), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 1,vjust=1, colour = "black"),
# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# guides(fill='none')
# ggsave(p,filename=paste0(save_dir,"/DBDs_affected.png"),width=7, height=4, dpi=400)



pdatagh <- data.frame(matrix(ncol=3, nrow=0))
for(k in 1:length(all_cancer)){

    tempdbd <- unique(dbd_purt_dt[dbd_purt_dt$CANCER == all_cancer[k],][,c(3,5)])
    wh <- which(tempdbd$ED %like% ',')
    if(length(wh) != 0){
        pdatam <- tempdbd[wh,]
        pdatau <- tempdbd[-wh,]
        for(i in 1:length(pdatam[[1]])){
            tempv <- unlist(strsplit(pdatam$ED[i], '[,]'))
            for(j in 1:length(tempv)){
                pdatau <- rbind(pdatau, data.frame(AS=pdatam$AS[i], ED=tempv[j]))
            }
        }
    }else{
        pdatau <- tempdbd
    }
    pdataux <- plyr::count(pdatau$ED)
    pdataux$CANCER <- rep(all_cancer[k], length(pdataux[[1]]))
    pdatagh <- rbind(pdatagh, pdataux)
}


# ##make the group "others (29)"
# wh <- which(pdatagh$x %in% pdata1$x)
# whe <- setdiff(seq(1,length(pdatagh[[1]])), wh)
# pdatagh$x[whe] <- paste0('Others (', length(pdata2[[1]]), ')')


##-- add number of unique splicing events ---
nnam <- c()
for(k in 1:length(pdatagh[[1]])){
    wh <- which(pdatauniy$x == pdatagh$x[k])
    nnam <- c(nnam, paste0(pdatauniy$x[wh],':',pdatauniy$freq[wh]))
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
pdatat <- data.frame(CANCER=tcancer, ED=tdbd, FRAC=tfrac)
pdatat$bval <- ifelse(pdatat$FRAC < 10, pdatat$FRAC, '>=10')
# pdatat$DBD <- factor(pdatat$DBD, levels=pdatax$x)


cols <- c('#800026','#000000','#ffffff','#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026')
# c('#08306b', '#6baed6','#deebf7', '#ffffcc','#fed976','#fd8d3c','#e31a1c','#800026')#'#2171b5', 
# cols <- rev(brewer.pal(11,"Spectral"))
basesize <- 8
p <- ggplot(pdatat, aes(ED, CANCER)) + geom_tile(aes(fill = bval), color='black')+scale_fill_manual(values=cols)
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
ggsave(p,filename=paste0(save_dir,'/Perturbed_DBDs_tile.png'),width=3.5, height=4, dpi=500)


##--- hypergeometruc test for each ED type ---------

bgsize <- sum(pdata$freq)
temps <- pdataq[pdataq$freq > 0, ]#alldata[alldata$FLAG == 'PTSEs', ]
sampleSize <- sum(temps$freq)
all_dbds <- unique(temps$x)
hyp1 <- c()
for(k in 1:length(all_dbds)){

    setB <- temps[temps$x == all_dbds[k], ]$freq
    setA <- pdata[pdata$x == all_dbds[k], ]$freq
    hyp1 <- c(hyp1, phyper(setB-1,setA,bgsize-setA, sampleSize))

}

hyp2 <- c()
for(k in 1:length(all_dbds)){

    setB <- temps[temps$x == all_dbds[k], ]$freq
    setA <- pdata[pdata$x == all_dbds[k], ]$freq
    hyp2 <- c(hyp2, phyper(setB-1,setA,bgsize-setA, sampleSize, lower.tail=FALSE))

}
hyp <- c(hyp1, hyp2)
hypx <- p.adjust(hyp, 'fdr')
##---------------------------------------------------


##---- alternative visualization of the fraction of perturbed DBD types ----
event_dbd_unq <- unique(dbd_purt_dt[, c(2,5)])
wh <- which(event_dbd_unq$ED %like% ',')
pdatamul <- event_dbd_unq[wh,]
pdatauni <- event_dbd_unq[-wh,]
for(k in 1:length(pdatamul[[1]])){
    tempv <- unlist(strsplit(pdatamul$ED[k], '[,]'))
    for(j in 1:length(tempv)){
        pdatauni <- rbind(pdatauni, data.frame(SYMBOL=pdatamul$SYMBOL[k], ED=tempv[j]))
    }
}
pdatauni <- unique(pdatauni)
pdataunix <- plyr::count(pdatauni$ED)
pdatauniy <- pdataunix
# wh <- which(pdataunix$x %in% pdata1$x)
# whe <- setdiff(seq(1,length(pdataunix[[1]])), wh)
# pdataunix$x[whe] <- paste0('Others (', length(pdata2[[1]]), ')')
# pdatauniy <- aggregate(pdataunix$freq, by=list(Category=pdataunix$x), FUN=sum)

pdatax1 <- pdata
pdatauniy1 <- pdatauniy
pdatax1$frac <- signif((pdatax1$freq/sum(pdatax1$freq)*100),3)
pdatauniy1$frac <- signif((pdatauniy1$freq/sum(pdatauniy1$freq)*100),3)

colnames(pdatauniy1) <- c('x','freq','frac')
pdatauniy1$FLAG <- rep('PTSEs',length(pdatauniy1[[1]]))
pdatax1$FLAG <- rep('Background',length(pdatax1[[1]]))
alldata <- rbind(pdatax1, pdatauniy1)

cols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
    '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#000000',
    '#ffffff','#636363','#bdbdbd')
p <- ggplot(alldata, aes(FLAG, frac, fill=x)) + 
geom_bar(stat="identity",position='stack', color='black')+
theme(legend.text=element_text(size=12))
basesize <- 10
p <- p + theme_bw() +
scale_x_discrete(name="") + 
scale_y_continuous(name="% of occurence/perturbation") +
scale_fill_manual(values=cols)+
# geom_text(aes(label=x), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
  theme(axis.text.x = element_text(size = basesize * 1,angle = 60, hjust = 1,vjust=1, colour = "black"), 
                         axis.text.y = element_text(size = basesize * 1, colour='black'), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=basesize*1.2),axis.title.y = element_text(size=basesize*1.2), 
        plot.title = element_text(size=basesize*1), strip.text.x = element_text(size = basesize * 1, colour = "black", angle = 0), 
        strip.text.y = element_text(size = basesize * 1, colour = "black", angle = 0), 
        legend.text=element_text(size=basesize), legend.title=element_text(size=basesize*1.2))+
guides(fill=guide_legend(title="Effector domain type"),ncol=2)
ggsave(p,filename=paste0(save_dir,"/EDs_affected_frac.png"),width=3.5, height=5, dpi=500)


#####------ overlap between TFs with DBD perturbations and TFs with ED perturbations-----
perturbed_tfs <- '../results_new/TF_splicing/PTSEs.xlsx'
perturbed_tfs_DBD <- '../results_new/TF_DBD/PTSEs_perturbing_DBDs.xlsx'
perturbed_tfs_ED <- '../results_new/TF_ED/PTSEs_perturbing_EDs.xlsx'
temp_dbd <- c()
temp_ed <- c()
temp_all <- c()
for(k in 1:length(all_cancer)){
    tdbd <- openxlsx::read.xlsx(perturbed_tfs_DBD, k)$AS_ID
    ted <- openxlsx::read.xlsx(perturbed_tfs_ED, k)$AS_ID
    if(k < 4){
        temps <- openxlsx::read.xlsx(perturbed_tfs, k)
    }else{
        temps <- openxlsx::read.xlsx(perturbed_tfs, k+1)
    }
    temp_all <- union(temp_all, temps$SYMBOL)
    temp_dbd <- union(temp_dbd, temps[temps$AS_ID %in% tdbd, ]$SYMBOL)
    temp_ed <- union(temp_ed, temps[temps$AS_ID %in% ted, ]$SYMBOL)
}


###---- plot for splicing events perturbing DBD, ED, both of DBD and ED, and none of DBD or ED -----------------

ed_dbd <- intersect(temp_dbd, temp_ed)
temp_oth <- setdiff(temp_all, union(temp_ed, temp_dbd))
temp_dbdx <- setdiff(temp_dbd, temp_ed)
temp_edx <- setdiff(temp_ed, temp_dbd)

pdata <- data.frame(Type=c('DBD', 'ED', 'DBD and ED','Other'), 
    countx=c(length(temp_dbdx),length(temp_edx), length(ed_dbd),length(temp_oth))
    )

pdata$count <- (pdata$countx/sum(pdata$countx))*100
pdata$count <- pdata$countx
library(scales)
library(ggrepel)
library(tidyverse)

pdata$count <- signif(pdata$count, 2)
df2 <- pdata %>% 
  mutate(csum = rev(cumsum(rev(count))), 
         pos = count/2 + lead(csum, 1),
         pos = if_else(is.na(pos), count/2, pos))

p <- ggplot(pdata, aes(x = "" , y = count, fill = fct_inorder(Type))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(count)),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "# of PTSE \naffected TFs"))+theme_void() +theme(legend.position='bottom')

ggsave(p,filename=paste0(save_dir,"/Pi_chart.png"),width=4, height=4, dpi=400)

#########--------------------------------------------------------------------------------


# ##--- overlap of the non DBD or ED perturbed TF with other TFs of PTSEs ----
# input_dir <- '../data/PSI_data'
# fdr <- 0.05
# all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
# not_included <- list()
# for(k in 1:length(all_cancer)){

#     temp <- data.table::fread(all_files[k], sep='\t')
#     whx <- which(toupper(temp$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
#     temp1 <- temp[whx,]

#     wha <- which(temp$FDR < fdr)
#     whb <- which(abs(temp$MEAN_NORMAL-temp$MEAN_CANCER) > fdr)
#     wh <- intersect(wha, whb)
#     tempx <- temp[wh, ]
#     whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
#     tempy <- tempx[whx,]
#     tempz <- tempy[tempy$symbol %in% setdiff(tempy$symbol, temp_oth),]
#     not_included[[k]] <- tempz[order(-abs(tempz$MEAN_DIFF)), ]

# }



