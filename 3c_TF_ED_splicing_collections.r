##############################################################################################
# Purpose: Add ED (effector domain) information
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)
options(warn=2)

input_dir <- '../data/PSI_data'
input_dirx <- '../data/uniprot_Ensembl_Exon_map_DBD_AS'
output_dir <- '../data/uniprot_Ensembl_Exon_map_DBD_ED_AS'

if(dir.exists(output_dir)){
    unlink(output_dir, recursive=TRUE)
}
dir.create(output_dir, recursive=TRUE)

all_filesx <- list.files(input_dirx, pattern='*.txt', full.names=TRUE)
alluni <- unlist(lapply(strsplit(basename(all_filesx), '[_]'),'[[',1))
all_uniprot <- unique(alluni)

ed_file1 <- as.data.frame(readxl::read_excel('../data/1-s2.0-S1097276521009576-mmc8.xlsx', 7))
ed_file1 <- ed_file1[, c(1,2)]
ed_file2 <- as.data.frame(readxl::read_excel('../data/1-s2.0-S1097276521009576-mmc8.xlsx', 5))
ed_file2 <- ed_file2[, c(1,2,3,4,5,6)]
## From this study: Compendium of human transcription factor effector domains, Molecular Cell 2022
ed_file <- merge(ed_file1, ed_file2, by='Effector domain ID')

tuniprot <- intersect(all_uniprot, ed_file$`Uniprot ID`)
## Of the 1,459 TFs with a matching Ensembl transcript, 536 have ED information ----

for(k in 1:length(all_uniprot)){

    tempuni <- all_filesx[which(alluni == all_uniprot[k])]
    temped <- ed_file[ed_file$`Uniprot ID` == all_uniprot[k], ]

    if(nrow(temped) == 0){

        for(j in 1:length(tempuni)){
            tempf <- data.table::fread(tempuni[j])
            tempf$ED <- rep('-',length(tempf[[1]]))
            data.table::fwrite(tempf,paste0(output_dir,'/',basename(tempuni[j])), sep='\t', row.names=FALSE, quote=FALSE)
        }

    }else{
        for(j in 1:length(tempuni)){
            tempf <- data.table::fread(tempuni[j])
            tempf$ED <- rep('-',length(tempf[[1]]))
            # tempf$ED_TF_Family <- rep('-',length(tempf[[1]]))
            # tempf$ED_Cluster_ID <- rep('-',length(tempf[[1]]))

            for(i in 1:length(temped[[1]])){
                tempedx <- temped[i,]
                start <- strsplit(tempedx$Coordinates, '[-]')[[1]][1]
                end <- strsplit(tempedx$Coordinates, '[-]')[[1]][2]
                tseq <- seq(start, end)
                wh <- which(tempf$UNIPROT_SEQ_NUM %in% tseq)
                tempf$ED[wh] <- tempedx$`Domain type`
                # tempf$ED_TF_Family[wh] <- tempedx$`TF Family`
                # tempf$ED_Cluster_ID[wh] <- tempedx$`Cluster ID`
            }
            
            data.table::fwrite(tempf,paste0(output_dir,'/',basename(tempuni[j])), sep='\t', row.names=FALSE, quote=FALSE)
        }
        
    }

    cat('Protein',k,'of',length(all_uniprot),'done\n')
}

