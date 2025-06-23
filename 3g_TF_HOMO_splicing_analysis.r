##############################################################################################
# Purpose: For each cancer type, see whether homodimer residues are removed by AS
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)
library(dorothea)
library(pheatmap)

save_dir <- '../results_new/TF_HOMODIMER'
if(dir.exists(save_dir)){
    unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

psi_input <- '../data/PSI_data'
as_input <- '../data/uniprot_Ensembl_Exon_map_DBD_ED_AS_HD'

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
            wht <- which(temp1$HOMODIMER_DEGREE != '-')
            if(length(wht) >= atleast_DBD){
                temp1 <- temp1[wht,]
                all_filesx1 <- c(all_filesx1, all_files[j])
                events1 <- union(events1, setdiff(unique(temp1$ES), '-'))
            }
        }

        if(nrow(temp2) > 0){
            wht <- which(temp2$HOMODIMER_DEGREE != '-')
            if(length(wht) >= atleast_DBD){
                temp2 <- temp2[wht,]
                all_filesx2 <- c(all_filesx2, all_files[j])
                events2 <- union(events2, setdiff(unique(temp2$AP), '-'))
            }
        }

        if(nrow(temp3) > 0){
            wht <- which(temp3$HOMODIMER_DEGREE != '-')
            if(length(wht) >= atleast_DBD){
                temp3 <- temp3[wht,]
                all_filesx3 <- c(all_filesx3, all_files[j])
                events3 <- union(events3, setdiff(unique(temp3$AT), '-'))
            }
        }

        if(nrow(temp4) > 0){
            wht <- which(temp4$HOMODIMER_DEGREE != '-')
            if(length(wht) >= atleast_DBD){
                temp4 <- temp4[wht,]
                all_filesx4 <- c(all_filesx4, all_files[j])
                events4 <- union(events4, setdiff(unique(temp4$AD), '-'))
            }
        }

        if(nrow(temp5) > 0){
            wht <- which(temp5$HOMODIMER_DEGREE != '-')
            if(length(wht) >= atleast_DBD){
                temp5 <- temp5[wht,]
                all_filesx5 <- c(all_filesx5, all_files[j])
                events5 <- union(events5, setdiff(unique(temp5$AA), '-'))
            }
        }

        if(nrow(temp6) > 0){
            wht <- which(temp6$HOMODIMER_DEGREE != '-')
            if(length(wht) >= atleast_DBD){
                temp6 <- temp6[wht,]
                all_filesx6 <- c(all_filesx6, all_files[j])
                events6 <- union(events6, setdiff(unique(temp6$ME), '-'))
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

##----------------------------------------------------------------------------


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


##--------------------------------------------------------------------

##--- Save all DBD perturbing events -----------------------------

##--- Which of the event pairs are also have affected DBDs -----
store_val <- function(y1){
    if(nrow(y1) == 0){
        flag <- 0
    }else{
        y2 <- y1[y1$HOMODIMER_DEGREE != '-', ]
        if(nrow(y2) < atleast_DBD){ ## consider an event to perturb a DBD if the event affects at least 2 DBD positions
            flag <- 0
        }else{
            flag <- unique(y2$HOMODIMER_DEGREE)
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
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Events_perturbing_HDs.xlsx'))

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
            } ## some gene symbols of TFs were not mapped to Ensembl
            tmap <- as.data.frame(data.table::fread(paste0(as_input,'/',tuni,'_',all_cancer[k],'.txt')))
            
            if(tempyy$splice_type[j] == 'AP'){
                y1 <- tmap[tmap$AP == tempyy$as_id[j], ]
                odbd <- c(odbd, store_val(y1)) 
                odbdt <- c(odbdt, store_val(y1)) 
            }else if(tempyy$splice_type[j] == 'AT'){
                y1 <- tmap[tmap$AT == tempyy$as_id[j], ]
                odbd <- c(odbd, store_val(y1)) 
                odbdt <- c(odbdt, store_val(y1))  
            }else if(tempyy$splice_type[j] == 'ES'){
                y1 <- tmap[tmap$ES == tempyy$as_id[j], ]
                odbd <- c(odbd, store_val(y1)) 
                odbdt <- c(odbdt, store_val(y1))   
            }else if(tempyy$splice_type[j] == 'AD'){
                y1 <- tmap[tmap$AD == tempyy$as_id[j], ]
                odbd <- c(odbd, store_val(y1)) 
                odbdt <- c(odbdt, store_val(y1))  
            }else if(tempyy$splice_type[j] == 'AA'){
                y1 <- tmap[tmap$AA == tempyy$as_id[j], ]
                odbd <- c(odbd, store_val(y1))  
                odbdt <- c(odbdt, store_val(y1)) 
            }else if(tempyy$splice_type[j] == 'ME'){
                y1 <- tmap[tmap$ME == tempyy$as_id[j], ]
                odbd <- c(odbd, store_val(y1)) 
                odbdt <- c(odbdt, store_val(y1)) 
            }else{
                odbd <- c(odbd, 0) 
                odbdt <- c(odbdt, 0) 
            }
            
            asid <- c(asid, tempyy$as_id[j])
            tcancer <- c(tcancer, all_cancer[k])
            corr1 <- c(corr1, tempyy$MEDIAN_CANCER[j])
            corr2 <- c(corr2, tempyy$MEDIAN_NORMAL[j])
            corr3 <- c(corr3, tempyy$MEDIAN_DIFF[j])
            corr4 <- c(corr4, tempyy$FDR[j])

            asidt <- c(asidt, tempyy$as_id[j])
            corr1t <- c(corr1t, tempyy$MEDIAN_CANCER[j])
            corr2t <- c(corr2t, tempyy$MEDIAN_NORMAL[j])
            corr3t <- c(corr3t, tempyy$MEDIAN_DIFF[j])
            corr4t <- c(corr4t, tempyy$FDR[j])
            tgene <- c(tgene, tempyy$symbol[j])
        }

    }
    
    ## save excel sheet ----
    tdatat <- data.frame(GENE=tgene, AS=asidt, DBD=odbdt, MEDIAN_CANCER=corr1t, MEDIAN_NORMAL=corr2t, MEDIAN_DIFF=corr3t, FDR=corr4t)
    tdatat <- tdatat[tdatat$DBD != 0, ]
    tdatat <- tdatat[order(abs(tdatat$MEDIAN_DIFF), decreasing=TRUE), ]
    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], tdatat)
    openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Events_perturbing_HDs.xlsx'), overwrite = T)
}

tdata <- data.frame(CANCER=tcancer, AS=asid, HOMODIMER_DEGREE=odbd, MEDIAN_CANCER=corr1, MEDIAN_NORMAL=corr2, MEDIAN_DIFF=corr3, FDR=corr4)
data.table::fwrite(tdata, paste0('../data/Events_perturbing_HD.txt'), row.names=FALSE, quote=FALSE, sep='\t')

##----------------------------------------------------------------
dbd_purt_dt <- tdata

##--- Events plots ------
Event_count <- c()
num_sig_events <- c()
sig_events <- list()
as_tags <- c('ES','AP','AT','AD','AA','ME')
AA <- c()
AD <- c()
AP <- c()
AT <- c()
ES <- c()
ME <- c()

for(k in 1:length(all_cancer)){

    tempf <- data.table::fread(all_filesxx[k], sep='\t')
    # tag <- c()
    ect <- c()
    for(j in 1:length(Events_list)){
        ect <- union(ect, Events_list[[j]][[k]])
        # tag <- c(tag, rep(as_tags[j], length(Events_list[[j]][[k]])))
    }
    Event_count <- c(Event_count, length(ect))

    tempx <- tempf[tempf$as_id %in% ect, ]

    num_sig_events <- c(num_sig_events, length(tempx[[1]]))

    sig_events[[k]] <- tempx$as_id

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
pdata <- data.frame(cancer=all_cancer, count=num_sig_events)
pdata_numbers <- pdata
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
maxv <- max(pdata$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of AS events perturbing \na DBD region of TFs", limits=c(0,maxv+20)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Events_perturbing_HDs.png"),width=3.5, height=3, dpi=400)


##--- plot the number of splicing events of different types affecting TFs -----------------
pData <- data.frame(A=AA, B=AD, C=AP, D=AT, E=ES, F=ME, X=num_sig_events, Cancer=all_cancer)
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
ggsave(ppx,filename=paste0(save_dir, "/Events_types_perturbing_HDs.png"),width=7, height=3, dpi=400)


###------ Distribution of perturbation values ----------------------------------------------------
##------------------------------------------------
tpos <- c()
tneg <- c()
tcancer <- c()
tdiff <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_filesxx[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    tempy <- tempx[tempx$as_id %in% all_events[[k]], ]

    tempy$POSP <- tempy$POS/paired_sam[[2]][k]
    tempy$NEGP <- tempy$NEG/paired_sam[[2]][k]

    tpos <- c(tpos, tempy$POSP)
    tneg <- c(tneg, tempy$NEGP)
    tcancer <- c(tcancer, rep(all_cancer[k],length(tempy[[1]])))
    tdiff <- c(tdiff, tempy$MEDIAN_DIFF)

}

pdata <- data.frame(CANCER=tcancer, High=tpos, Low=tneg, impact=tdiff)
##--- choose the max of # of patients in which an event is gained or lost ----
flag <- c()
maxp <- c()
for(k in 1:length(pdata[[1]])){
    temp <- pdata[k,]
    maxp <- c(maxp, max(temp$High, temp$Low))
    if(temp$High > temp$Low){
        flag <- c(flag, 'High')
    }else{
        flag <- c(flag, 'Low')
    }
}

pdata$MAXP <- maxp
pdata$FLAG <- flag

p <- ggplot(pdata, aes(MAXP, impact, color=FLAG)) + 
geom_point(size=0.4)+geom_density_2d(colour='white', size=0.2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_y_continuous(name="Difference between \nmedian PSIs of paired samples", limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) + 
scale_x_continuous(name="Fraction of patients in a cancer type", limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
# geom_hline(yintercept=0.2, color='black', linetype='dashed', size=0.2)+
# geom_hline(yintercept=-0.2, color='black', linetype='dashed', size=0.2)+
scale_color_manual(values=c('#e41a1c','#377eb8'))+
guides(color=guide_legend(title="PSI in cancer\nvs.\npaired normal\nsamples",ncol=1, override.aes = list(size=2)))
ggsave(p,filename=paste0(save_dir,"/Perturbed_events_patients_HD.png"),width=4.5, height=3, dpi=600)


##-----------------------------------------------------------------------------------------

##-- Look into the splicing behaviour of the DBD perturbing events -------------------
cancer <- c()
mcr <- c()
mnl <- c()
asid <- c()
gene <- c()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_filesxx[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    tempy <- tempx[tempx$as_id %in% all_events[[k]], ]

    mcr <- c(mcr, tempy$MEDIAN_CANCER)
    mnl <- c(mnl, tempy$MEDIAN_NORMAL)
    cancer <- c(cancer, rep(all_cancer[k], length(tempy[[1]])))
    asid <- c(asid, tempy$as_id)
    gene <- c(gene, tempy$symbol)
}

pdata <- data.frame(cancer=cancer, gene=gene, asd=asid, median_cancer=mcr, median_normal=mnl)
pdata$median_diff <- pdata$median_cancer-pdata$median_normal

p <- ggplot(pdata, aes(cancer, median_diff)) + 
geom_jitter(aes(color=cancer),size=0.4)+
geom_violin()+
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
ggsave(p,filename=paste0(save_dir,"/Sig_events_TFs_HD_dist.png"),width=5, height=3, dpi=400)




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
        temp_ovl <- sig_events[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- intersect(temp_ovl, sig_events[[which(all_cancer == temp_comb[[i]][j])]])
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
maxv <- max(pdata$count)
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of splicing events", limits=c(0,maxv+20)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Events_perturbing_HDs_overlap.png"),width=3.5, height=3, dpi=400)



##--- plot the number of TFs -----------------
pdata <- data.frame(cancer=all_cancer, count=lengths(all_tfs))
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of TFs with DBDs \naffected by a splicing event", limits=c(0,max(pdata$count)+5)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/TFs_with_affected_HDs.png"),width=3.5, height=3, dpi=400)


##----- PSI profile of the pancancer (at least 10 cancer types) events -------------
events_to_consider <- num_events_combs[[10]]
tcancer <- c()
tvalue <- c()
tevent <- c()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_filesxx[k], sep='\t')
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
pdx <- t(as.matrix(pdata_mat))

color <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
breaks <- seq(-1, 1, 0.2)
p <- pheatmap(pdx,fontsize=3, color=color, breaks=breaks, cluster_rows=FALSE, cluster_cols=FALSE,cellheight=5, cellwidth = 5)
ggsave(p,filename=paste0(save_dir, "/Pancancer_events_HD.png"),width=7, height=5, dpi=600)


