##############################################################################################
# Purpose: Add ED (effector domain) information
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)
library(Biostrings)
library(ggrepel)
options(warn=2)

input_dirx <- 'data/uniprot_Ensembl_Exon_map_DBD_ED_AS'
output_dir <- 'data/uniprot_Ensembl_Exon_map_DBD_ED_AS_Lambourne'
save_dir <- 'results_rep/Lambourne'

if(dir.exists(save_dir)){
    unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)


if(dir.exists(output_dir)){
    unlink(output_dir, recursive=TRUE)
}
dir.create(output_dir, recursive=TRUE)

tf_ensemb_map <- as.data.frame(data.table::fread('data/TF_ensembl_uniprot.txt', sep='\t'))
tf_ensemb_map$Ensembl_protein_id <- unlist(lapply(strsplit(tf_ensemb_map$Ensembl_protein_id, '[.]'), '[[',1))

all_filesx <- list.files(input_dirx, pattern='*.txt', full.names=TRUE)
alluni <- unlist(lapply(strsplit(basename(all_filesx), '[_]'),'[[',1))
all_uniprot <- unique(alluni)
all_hgnc <- tf_ensemb_map$HGNC_symbol[which(tf_ensemb_map$Uniprotswissprot %in% all_uniprot)]

##-- Known isoform list with experimental functional information ----
##-- Taken from this study (Data S1): Widespread variation in molecular interactions and regulatory properties among transcription factor isoforms

lamf <- data.table::fread('data/SuppTable_CloneList.txt')
lamf_tfs <- unique(lamf$gene_symbol)

my_lamf_tfs <- intersect(lamf_tfs, all_hgnc)

my_uni_lamf <- tf_ensemb_map$Uniprotswissprot[which(tf_ensemb_map$HGNC_symbol %in% my_lamf_tfs)]

all_filesy <- all_filesx[which(alluni %in% my_uni_lamf)]
alluniy <- unlist(lapply(strsplit(basename(all_filesy), '[_]'),'[[',1))

counterxx <- 0
temp_isoforms <- c()
for(k in 1:length(all_filesy)){

    tempf <- data.table::fread(all_filesy[k])
    temp_hgnc <- tf_ensemb_map$HGNC_symbol[which(tf_ensemb_map$Uniprotswissprot == alluniy[k])]
    temp_lamb <- lamf[lamf$gene_symbol == temp_hgnc, ]
    toalign1 <- paste(tempf$UNIPROT, collapse='')

    lam_seq_names <- c()
    afpos <- vector(mode = "list", length = length(temp_lamb[[1]]))
    counterxx <- counterxx + nrow(temp_lamb)
    temp_isoforms <- union(temp_isoforms, temp_lamb$clone_id)

    for(j in 1:length(temp_lamb[[1]])){ ###-- for each isoform ----

        lam_seq_names <- c(lam_seq_names, temp_lamb$clone_id[j])

        toalign2 <- temp_lamb$aa_seq[j]
        xxx <- pwalign::pairwiseAlignment(pattern=toalign1, subject=toalign2,
             substitutionMatrix = "BLOSUM62", 
             gapOpening = -2,
             gapExtension = -8, 
             scoreOnly = FALSE)

        xx1 <- unlist(strsplit(as.character(pwalign::pattern(xxx)), split=''))
        xx2 <- unlist(strsplit(as.character(pwalign::subject(xxx)), split=''))

        temp_pos <- c()
        counter <- 0
        for(i in 1:length(xx1)){
            if(xx1[i] != '-' & xx2[i] != '-'){
                counter <- counter+1
            }else if(xx1[i] != '-' & xx2[i] == '-'){
                counter <- counter+1
                temp_pos <- c(temp_pos, counter)
            }
        }

        tempst <- rep(0,nchar(toalign1))
        if(length(temp_pos) != 0){tempst[temp_pos] <- 1}
        afpos[[j]] <- tempst

    }

    afposx <- as.data.frame(afpos)
    colnames(afposx) <- lam_seq_names

    tempg <- cbind(tempf, afposx)
    data.table::fwrite(tempg,paste0(output_dir,'/',basename(all_filesy[k])), sep='\t', row.names=FALSE, quote=FALSE)

    cat('Map',k,'of',length(all_filesy),'done\n')

}


###--- analyze the data to see which isoforms could be impacted by which AS event---
### Do these findings agree with the experimental findings??
psi_input <- 'data/PSI_data'
all_filesxx <- gtools::mixedsort(list.files(psi_input, pattern='*PSI_paired.txt', full.names=TRUE))
all_filesxx <- all_filesxx[-4]
all_cancer <- substr(basename(all_filesxx), 1,4)
ovr_tfs <- c()
tcancer <- c()
tevent <- c()
tisoform <- c()
teid <- c()
not_c <- c()
not_e <- c()

for(k in 1:length(all_cancer)){

    all_files <- list.files(output_dir, pattern=paste0('*',all_cancer[k],'.txt'), full.names=TRUE)
    af_tfs <- 0

    for(j in 1:length(all_files)){

        temp <- data.table::fread(all_files[j])
        evpos <- vector(mode = "list", length = 6)

        evpos[[1]] <- which((temp$ES != '')&(!is.na(temp$ES)))
        evpos[[2]] <- which((temp$AP != '')&(!is.na(temp$AP)))
        evpos[[3]] <- which((temp$AT != '')&(!is.na(temp$AT)))
        evpos[[4]] <- which((temp$AD != '')&(!is.na(temp$AD)))
        evpos[[5]] <- which((temp$AA != '')&(!is.na(temp$AA)))
        evpos[[6]] <- which((temp$ME != '')&(!is.na(temp$ME)))

        wht <- which(lengths(evpos) > 0)

        if(length(wht) > 0){ ##-- if at least one event type is present
            af_tfs <- af_tfs+1 ## Number of TFs wiht at least one event

            for(i in 1:length(wht)){ ##-- for each event type --

                loop <- seq(22, length(temp))
                for(h in 1:length(loop)){ ##-- for each isoform
                    tempx1 <- temp[[loop[h]]][evpos[[wht[i]]]] ## positions of the uniprot protein missing (1) from the isoform and overlapping with the event
                    tempx2 <- temp[[loop[h]]][setdiff(seq(1,length(temp[[1]])),evpos[[wht[i]]])]
                    whx1 <- which(tempx1 == 0)
                    whx2 <- which(tempx2 == 1)
                    if(length(whx1) == 0 & length(whx2) == 0){
                        tcancer <- c(tcancer, all_cancer[k])
                        tevent <- c(tevent, paste(unique(temp[[14+wht[i]]][evpos[[wht[i]]]]), collapse=';') )
                        tisoform <- c(tisoform, colnames(temp)[loop[h]])
                        teid <- c(teid, colnames(temp)[14+wht[i]])
                    }else{
                        not_c <- c(not_c, k)
                        not_e <- c(not_e, j)
                    }

                }
                # temp_pos <- evpos[[wht[i]]] 
            }
        }
    }
    ovr_tfs <- c(ovr_tfs, af_tfs) ## number of TFs among the 223 proteins mapped to Lambourne, where at least one significant event is mapped...

}

pdata <- data.frame(CANCER=tcancer, AS_ID=tevent, ISOFORM=tisoform, SPLICE_TYPE=teid)
wh1 <- which(pdata$AS_ID %like% ';')
wh2 <- setdiff(seq(1, length(pdata[[1]])), wh1)
pdata1 <- pdata[wh2,]
tcancer <- c()
tevent <- c()
tisoform <- c()
teid <- c()
for(k in 1:length(wh1)){
    tids <- unlist(strsplit(as.character(pdata[wh1[k], ]$AS_ID), '[;]'))
    tdata <- data.frame(CANCER=rep(pdata[wh1[k], ]$CANCER, length(tids)), AS_ID=tids, 
        ISOFORM=rep(pdata[wh1[k], ]$ISOFORM, length(tids)), SPLICE_TYPE=rep(pdata[wh1[k], ]$SPLICE_TYPE, length(tids)))
    pdata1 <- rbind(pdata1, tdata)
}


##--- overlap with the transcriptional activity ----------------
temp_genes <- unique(unlist(lapply(strsplit(pdata$ISOFORM,'-'),'[[', 1)))
tempf <- data.table::fread('data/Data_S3.tsv')
pdatax <- tempf[tempf$gene_symbol %in% temp_genes, ]
activity_logfc <- c()
for(k in 1:length(pdata1[[1]])){

    tempg <- strsplit(pdata1$ISOFORM[k],'[-]')[[1]][1]
    tempx <- lamf[lamf$gene_symbol == tempg, ][,c(1,2,3)]
    whx <- which(tempx$isoform_status %like% 'reference')
    ref_iso <- tempx$clone_id[whx]
    why <- which(pdatax$clone_id == ref_iso)
    whz <- which(pdatax$clone_id == as.character(pdata1$ISOFORM[k]))

    if((length(why) > 0) & (length(whz) > 0)){
        ## determine fold change --
        activity_logfc <- c(activity_logfc, (pdatax[why, ]$M1H_mean)/(pdatax[whz, ]$M1H_mean))

    }else{
        activity_logfc <- c(activity_logfc, NA)
    }

}
pdata1$ACTVITY <- activity_logfc

tsheets <- readxl::excel_sheets('results_rep/TF_splicing/Supplementary_Table_S2.xlsx')
dpsi <- c()
for(k in 1:length(pdata1[[1]])){

    indx <- which(tsheets == pdata1$CANCER[k])
    temp_splice <-  as.data.frame(readxl::read_excel('results_rep/TF_splicing/Supplementary_Table_S2.xlsx',indx))
    dpsi <- c(dpsi, temp_splice[temp_splice$AS_ID == pdata1$AS_ID[k], ]$MEAN_DIFF)

}

pdata1$MEAN_DIFF <- dpsi


##---- plot delta psi vs transcription activity fold change ----
pdata2 <- pdata1[complete.cases(pdata1), ]

##--- remove clone-ids marked as reference in lambourne et al. ----
wh <- which(lamf$isoform_status %like% "reference")
whx <- which(pdata2$ISOFORM %in% lamf$clone_id[wh])
pdata2 <- pdata2[-whx, ]
##-----------------------------------------------------------------

pdata2$SYMBOL <- unlist(lapply(strsplit(pdata2$ISOFORM, '[-]'), '[[', 1))
pdata2$ID <- paste0(pdata2$SYMBOL,'_',pdata2$AS_ID,'_',pdata2$SPLICE_TYPE,'\n(',pdata2$CANCER,',',' \u0394PSI: ', signif(pdata2$MEAN_DIFF,3),')')
pdata2$PL <- ''
# wh1 <- which(pdata2$AS_ID %in% c('76094','77374'))
# wh2 <- which(abs(pdata2$MEAN_DIFF) > 0.25)
# wh <- intersect(wh1, wh2)
# pdata2$PL[wh] <- pdata2$ID[wh]

# whx1 <- which(pdata2$CANCER %in% c('UCEC','KIRC'))
# whx2 <- which(pdata2$SYMBOL %in% c('HSF2', 'NFYA'))
# whx <- intersect(whx1, whx2)
whx <- c(8,24,23,16,22)
pdata2$PL[whx] <- pdata2$ID[whx]


p <- ggplot(pdata2, aes(ACTVITY, MEAN_DIFF, color=SPLICE_TYPE, label=PL)) + 
# p <- ggplot(pdata2, aes(ACTVITY, MEAN_DIFF, color=SPLICE_TYPE)) + 
geom_point()+
theme(legend.text=element_text(size=12))
basesize <- 8
p <- p + 
scale_x_continuous(name="Transcription activity fold change (alt/ref)", limits=c(-1.5,2.5), breaks = seq(-1.5,2.5, by = 0.5)) + 
scale_y_continuous(name="Mean \u0394PSI", limits=c(-0.6, 0.25)) +
scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d'))+
# geom_hline(yintercept=0, linetype='dashed')+
# geom_vline(xintercept=1, linetype='dashed')+
geom_text_repel(family = "Poppins",
    max.overlaps=Inf,
                      size = 2.6,
                      color='black',
                      arrow = arrow(length = unit(0.010, "npc")),
                      box.padding=0.5,
                      min.segment.length = 0) +
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),axis.title=element_text(size=1*basesize),
    legend.text=element_text(size=basesize), legend.title=element_text(size=basesize),
    panel.border = element_blank())+
guides(color=guide_legend(title="Type of\nalternative\nsplicing\nevent", ncol=1))
ggsave(p,filename=paste0(save_dir,"/trancriptional_activity.png"),width=3.5, height=3, dpi=500)


##--- save excel file ------
pdata3 <- pdata2[,c(1,7,2,4,6,3,5)]
colnames(pdata3) <- c('CANCER','SYMBOL','AS_ID','AS_TYPE','MEAN_DIFF','ISOFORM','ACTIVITY_CHANGE')
wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Supplementary_Table_S9.xlsx'))
##--- column name explanations --
tdatat <- data.frame(`Table fields explanation`=c('Supplementary Table S9: List of PTSEs leading to TF isoforms with altered transcriptional activity as per Lambourne et al.',
    'Explanation of the table fields are outlined below','CANCER: TGCA cancer code',
    'SYMBOL: HGNC gene symbol',
   'AS_ID: Alternative splicing ID provided by TCGA SpliceSeq',
   'AS_TYPE: Type of alternative splicing',
   'MEAN_DIFF: Mean of the differences of PSI values between paired normal and cancer samples',
   'ISOFORM: Isoform from the Lambourne at al. study that overlapped with this study',
   'ACTIVITY_CHANGE: M1H activity fold chnage between ISOFORM and the corresponding reference protein taken from Lambourne et al.'
     ))
openxlsx::addWorksheet(wb1, sheetName = 'INDEX')
openxlsx::writeData(wb1, sheet = 'INDEX', tdatat, colNames=FALSE)
openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S9.xlsx'), overwrite = T)

openxlsx::addWorksheet(wb1, sheetName = 'Table S9')
openxlsx::writeData(wb1, sheet = 'Table S9', pdata3)
openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S9.xlsx'), overwrite = T)



# ##--- determine which splicing event to show details for -------------------------
# perturbed_tfs_unq <- 'results_rep/TF_splicing/Supplementary_Table_S3.xlsx'
# perturbed_tfs_DBD <- 'results_rep/TF_DBD/Supplementary_Table_S4.xlsx'
# perturbed_tfs_ED <- 'results_rep/TF_ED/Supplementary_Table_S5.xlsx'
# perturbed_tfs_lam <- 'results_rep/Lambourne/Supplementary_Table_S8.xlsx'
# temp_dbd <- list()
# temp_ed <- list()
# loop <- length(all_cancer)+1
# for(k in 2:loop){
#     temp_dbd[[k-1]] <- openxlsx::read.xlsx(perturbed_tfs_DBD, k)$AS_ID
#     temp_ed[[k-1]] <- openxlsx::read.xlsx(perturbed_tfs_ED, k)$AS_ID
# }

# temp_uq <- list()
# for(k in 2:loop){
#     temp <- openxlsx::read.xlsx(perturbed_tfs_unq, k)
#     temp_uq[[k-1]] <- temp[!is.na(temp$CI_Type),]$AS_ID
# }

# temp_lam <- openxlsx::read.xlsx(perturbed_tfs_lam, 2)

# # temp_lam <- list()
# # for(k in 2:loop){
# #     temp_lam[[k-1]] <- temp[!is.na(temp$CI_Type),]$AS_ID
# # }

# ed_ovr <- list()
# dbd_ovr <- list()
# for(k in 1:length(all_cancer)){
#     tempc <- temp_lam[temp_lam$CANCER == all_cancer[k], ]

#     if(nrow(tempc) != 0){
#         ed_ovr[[k]] <- intersect(temp_ed[[k]], tempc$AS_ID)
#         dbd_ovr[[k]] <- intersect(temp_dbd[[k]], tempc$AS_ID)
#     }
# }




# # ##--- overlap with the protein-DNA interaction --------------------
# # tempf <- data.table::fread('data/SuppTable_eY1HResults.txt')
# # for(k in 1:length(pdata1[[1]])){

# #     tempg <- strsplit(pdata1$ISOFORM[k],'[-]')[[1]][1]
# #     tempx <- lamf[lamf$gene_symbol == tempg, ][,c(1,2,3)]
# #     whx <- which(tempx$isoform_status %like% 'reference')
# #     ref_iso <- tempx$clone_id[whx]
# #     why <- which(tempf$clone_id == ref_iso)
# #     whz <- which(tempf$clone_id == as.character(pdata1$ISOFORM[k]))

# #     if((length(why) > 0) & (length(whz) > 0)){
        
# #         ## number of protein-DNA interactions for the reference
# #         tpdna <- which(unlist(tempf[why,]) == 'TRUE')
# #         apdna <- which(unlist(tempf[whz,]) == 'TRUE')

# #     }else{
        
# #     }

# # }









