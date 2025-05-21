##############################################################################################
# Purpose: Add ED (effector domain) information
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)

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

ed_file <- as.data.frame(readxl::read_excel('../data/1-s2.0-S1097276521009576-mmc8.xlsx', 5))
## From this study: Compendium of human transcription factor effector domains, Molecular Cell 2022

tuniprot <- intersect(all_uniprot, ed_file$`Uniprot ID`)
## Of the 1,460 TFs with DBD information, 536 have ED information ----

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
            start <- strsplit(temped$Coordinates, '[-]')[[1]][1]
            end <- strsplit(temped$Coordinates, '[-]')[[1]][2]
            tseq <- seq(start, end)
            wh <- which(tempf$UNIPROT_SEQ_NUM %in% tseq)
            tempf$ED[wh] <- temped$`Domain type`
            data.table::fwrite(tempf,paste0(output_dir,'/',basename(tempuni[j])), sep='\t', row.names=FALSE, quote=FALSE)
        }
        
    }

    cat('Protein',k,'of',length(all_uniprot),'done\n')
}

