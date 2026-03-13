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
# sk <- c()
# sj <- c()

for(k in 1:length(all_filesy)){

    tempf <- data.table::fread(all_filesy[k])
    temp_hgnc <- tf_ensemb_map$HGNC_symbol[which(tf_ensemb_map$Uniprotswissprot == alluniy[k])]
    temp_lamb <- lamf[lamf$gene_symbol == temp_hgnc, ]
    toalign1 <- paste(tempf$UNIPROT, collapse='')

    lam_seq_names <- c()
    afpos <- vector(mode = "list", length = length(temp_lamb[[1]]))

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

        # if(nchar(toalign1) < nchar(toalign2)){
        #     sk <- c(sk, k)  
        #     sj <- c(sj, j) 
        # }
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
            af_tfs <- af_tfs+1
            # print(j)
            # break

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
                temp_pos <- evpos[[wht[i]]]
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
        activity_logfc <- c(activity_logfc, pdatax[whz, ]$M1H_mean/pdatax[why, ]$M1H_mean)

    }else{
        activity_logfc <- c(activity_logfc, NA)
    }

}
pdata1$ACTVITY <- activity_logfc

tsheets <- readxl::excel_sheets('results_rep/TF_splicing/Supplementary_file_S1.xlsx')
dpsi <- c()
for(k in 1:length(pdata1[[1]])){

    indx <- which(tsheets == pdata1$CANCER[k])
    temp_splice <-  as.data.frame(readxl::read_excel('results_rep/TF_splicing/Supplementary_file_S1.xlsx',indx))
    dpsi <- c(dpsi, temp_splice[temp_splice$AS_ID == pdata1$AS_ID[k], ]$MEAN_DIFF)

}

pdata1$MEAN_DIFF <- dpsi


##---- plot delta psi vs transcription activity fold change ----
pdata2 <- pdata1[complete.cases(pdata1), ]

p <- ggplot(pdata2, aes(ACTVITY, MEAN_DIFF, color=SPLICE_TYPE)) + 
geom_point()+
theme(legend.text=element_text(size=12))
basesize <- 10
p <- p + 
scale_x_continuous(name="Transcription activity fold change (ref/alt)", limits=c(-1,2), breaks = seq(-1,2, by = 0.2)) + 
scale_y_continuous(name="Mean \u0394PSI", limits=c(-0.6, 0.25)) +
scale_color_manual(values=c('#e41a1c','#377eb8'))+
geom_hline(yintercept=0, linetype='dashed')+
geom_vline(xintercept=0, linetype='dashed')+
# geom_text_repel(family = "Poppins",
#     max.overlaps=Inf,
#                       size = 3,
#                       color='black',
#                       arrow = arrow(length = unit(0.010, "npc")),
#                       min.segment.length = 0) +
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),axis.title=element_text(size=1*basesize),
    legend.text=element_text(size=basesize), legend.title=element_text(size=basesize),
    panel.border = element_blank())+
guides(color=guide_legend(title="Type of\nalternative\nsplicing\nevent", ncol=1))
ggsave(p,filename=paste0(save_dir,"/trancriptional_activity.png"),width=5, height=3, dpi=500)



# ##--- overlap with the protein-DNA interaction --------------------
# tempf <- data.table::fread('data/SuppTable_eY1HResults.txt')
# for(k in 1:length(pdata1[[1]])){

#     tempg <- strsplit(pdata1$ISOFORM[k],'[-]')[[1]][1]
#     tempx <- lamf[lamf$gene_symbol == tempg, ][,c(1,2,3)]
#     whx <- which(tempx$isoform_status %like% 'reference')
#     ref_iso <- tempx$clone_id[whx]
#     why <- which(tempf$clone_id == ref_iso)
#     whz <- which(tempf$clone_id == as.character(pdata1$ISOFORM[k]))

#     if((length(why) > 0) & (length(whz) > 0)){
        
#         print(k)

#     }else{
        
#     }

# }









