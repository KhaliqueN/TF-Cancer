##############################################################################################
# Purpose: For each cancer type, see whether DBDs are removed by AS
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)

save_dir <- '../results'
input_dir <- '../data/PSI_data'
input_dirx <- '../data/uniprot_Ensembl_Exon_mapx'
output_dir <- '../data/uniprot_Ensembl_Exon_map_DBD_AS'
# if(dir.exists(output_dir)){
#     unlink(output_dir, recursive=TRUE)
# }
# dir.create(output_dir)

fdr <- 0.05
diff <- 0.1
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
tf_ensemb_map <- as.data.frame(data.table::fread('../data/TF_ensembl_uniprot.txt', sep='\t'))

tcga_map <- data.table::fread(paste0(input_dir,'/TCGA_SpliceSeq_Gene_Structure.txt'))

all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI.txt', full.names=TRUE))
all_filesx <- list.files(input_dirx, pattern='*.txt', full.names=TRUE)
all_cancer <- substr(basename(all_files), 1,4)
## Some of the event coordinates frm TCGA splice seq does not overlap to any of the exons mapped from the canonical protein
## sequence from Uniprot. This is because the canonical seqeunce is not always the longest.


###--- EXON SKIPPING -----------------------------------------------------------------------------
for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    tempy <- tempx[which(toupper(tempx$symbol) %in% tfs$Gene_Symbol), ]
    tempz <- tempy[tempy$splice_type == 'ES', ]
    all_symbolx <- gtools::mixedsort(unique(tempz$symbol))

    ##--- only the processed TFs --- 1-to-1 mapped transcripts to TFs
    tid <- unlist(lapply(strsplit(basename(all_filesx), '[.]'),'[[',1))
    eid <- tf_ensemb_map[which(tf_ensemb_map$Uniprotswissprot %in% tid),]$Ensembl_gene_id
    back_symbols <- tfs[which(tfs$Ensembl_Gene_ID %in% eid), ]$Gene_Symbol 
    ## Out of 1639 curated TFs, 1460 are also present in the exon mapping file
    ### So, the splice event information will be stored for these 1460 TFs
    
    all_symbol <- intersect(all_symbolx, back_symbols)
    rest_symbols <- setdiff(back_symbols, all_symbol)

    for(j in 1:length(all_symbol)){ ## for k=1, you can check for j=6 and i=2

        tempz1 <- tempz[tempz$symbol == all_symbol[j], ]

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == all_symbol[j])])]

        temp_map <- data.table::fread(paste0(input_dirx,'/',temp_uniprot,'.txt'))

        nt_start <- temp_map$NT1 ## for each AA position
        nt_end <- temp_map$NT3
        ##------------------------------------------------------------------------
        as_event <- rep('-', length(temp_map[[1]]))
        temp_cancer <- rep('-', length(temp_map[[1]]))

        for(jj in 1:length(tempz1[[1]])){ ## for each splicing event

            tempz2 <- tempz1[jj, ]
            skipped_exons <- strsplit(tempz2$exons,'[:]')[[1]]

            for(i in 1:length(skipped_exons)){ ## for each skipped exon, store the ES event ID

                tcga_temp <- tcga_map[tcga_map$Symbol == all_symbol[j], ]
                tcga_temp <- tcga_temp[tcga_temp$Exon == skipped_exons[i], ]
                temp_strand <- tcga_temp$Strand

                for(h in 1:length(nt_start)){ ## for each dbd start
                    ## dbd start and end should be within tcga_temp start and end, in order to consider the dbd to be alternatively spliced

                    if(temp_strand == '-'){
                        wh1 <- (nt_start[h] >= tcga_temp$Chr_Stop)
                        wh2 <- (nt_end[h] <= tcga_temp$Chr_Start)
                    }else{
                        wh1 <- (nt_start[h] >= tcga_temp$Chr_Start)
                        wh2 <- (nt_end[h] <= tcga_temp$Chr_Stop)
                    }
                    
                    if(wh1 & wh2){ 
                        ##-- store the event id at the protein positions not included in the protein product of the sample with lower PSI----
                        as_event[h] <- tempz2$as_id
                        # print(paste0(j,':',i))
                    }
                }

            }
        }

        temp_map$ES <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    for(j in 1:length(rest_symbols)){

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == rest_symbols[j])])]

        temp_map <- data.table::fread(paste0(input_dirx,'/',temp_uniprot,'.txt'))
        as_event <- rep('-', length(temp_map[[1]]))
        temp_map$ES <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    cat('Cancer',k,'of',length(all_cancer),'done\n')
}


###--- ALTERNATIVE PROMOTER ---------------------------------------
for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    tempy <- tempx[which(toupper(tempx$symbol) %in% tfs$Gene_Symbol), ]
    tempz <- tempy[tempy$splice_type == 'AP', ]
    all_symbolx <- gtools::mixedsort(unique(tempz$symbol))

    ##--- only the processed TFs --- 1-to-1 mapped transcripts to TFs
    tid <- unlist(lapply(strsplit(unlist(lapply(strsplit(basename(all_filesx), '[_]'),'[[',1)), '[.]'),'[[',1))
    eid <- tf_ensemb_map[which(tf_ensemb_map$Uniprotswissprot %in% tid),]$Ensembl_gene_id
    back_symbols <- tfs[which(tfs$Ensembl_Gene_ID %in% eid), ]$Gene_Symbol 
    ## Out of 1639 curated TFs, 1460 are also present in the exon mapping file
    ### So, the splice event information will be stored for these 1460 TFs
    
    all_symbol <- intersect(all_symbolx, back_symbols)
    rest_symbols <- setdiff(back_symbols, all_symbol)

    for(j in 1:length(all_symbol)){ ## for k=1, you can check for j=6 and i=2

        tempz1 <- tempz[tempz$symbol == all_symbol[j], ]

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == all_symbol[j])])]

        temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))

        nt_start <- temp_map$NT1 ## for each AA position
        nt_end <- temp_map$NT3
        ##------------------------------------------------------------------------
        as_event <- rep('-', length(temp_map[[1]]))
        temp_cancer <- rep('-', length(temp_map[[1]]))

        for(jj in 1:length(tempz1[[1]])){ ## for each splicing event

            tempz2 <- tempz1[jj, ]
            skipped_exons <- tempz2$exons ## the coding protein starts from this exon
            ## meaning that the exons before this exon are missing from PSI sample
            ## mark positions before this exon as affected by this splice event

            for(i in 1:length(skipped_exons)){ ## for each skipped exon, store the ES event ID

                tcga_temp <- tcga_map[tcga_map$Symbol == all_symbol[j], ]
                tcga_temp <- tcga_temp[tcga_temp$Exon == skipped_exons[i], ]
                temp_strand <- tcga_temp$Strand

                for(h in 1:length(nt_start)){ ## for each dbd start
                    ## dbd start and end should be within tcga_temp start and end, in order to consider the dbd to be alternatively spliced

                    if(temp_strand == '-'){ ## if negative strand then genomic positions higher than the affect exon is not included
                        ## meaning anything higher than the start nt position of the affected exon
                        whx <- (nt_start[h] >= tcga_temp$Chr_Start)
                    }else{
                        whx <- (nt_start[h] <= tcga_temp$Chr_Start)
                    }
                    
                    if(whx){ 
                        ##-- store the event id at the protein positions not included in the protein product of the sample with higher PSI ----
                        as_event[h] <- tempz2$as_id
                    }
                }

            }
        }

        temp_map$AP <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    for(j in 1:length(rest_symbols)){

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == rest_symbols[j])])]

        temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))
        as_event <- rep('-', length(temp_map[[1]]))
        temp_map$AP <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    cat('Cancer',k,'of',length(all_cancer),'done\n')
}

###--- Work from here ---------------------------------------------

###--- ALTERNATIVE TERMINATOR ---------------------------------------
for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    tempy <- tempx[which(toupper(tempx$symbol) %in% tfs$Gene_Symbol), ]
    tempz <- tempy[tempy$splice_type == 'AT', ]
    all_symbolx <- gtools::mixedsort(unique(tempz$symbol))

    ##--- only the processed TFs --- 1-to-1 mapped transcripts to TFs
    tid <- unlist(lapply(strsplit(unlist(lapply(strsplit(basename(all_filesx), '[_]'),'[[',1)), '[.]'),'[[',1))
    eid <- tf_ensemb_map[which(tf_ensemb_map$Uniprotswissprot %in% tid),]$Ensembl_gene_id
    back_symbols <- tfs[which(tfs$Ensembl_Gene_ID %in% eid), ]$Gene_Symbol 
    ## Out of 1639 curated TFs, 1460 are also present in the exon mapping file
    ### So, the splice event information will be stored for these 1460 TFs
    
    all_symbol <- intersect(all_symbolx, back_symbols)
    rest_symbols <- setdiff(back_symbols, all_symbol)

    for(j in 1:length(all_symbol)){ ## for k=1, you can check for j=6 and i=2

        tempz1 <- tempz[tempz$symbol == all_symbol[j], ]

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == all_symbol[j])])]

        temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))

        nt_start <- temp_map$NT1 ## for each AA position
        nt_end <- temp_map$NT3
        ##------------------------------------------------------------------------
        as_event <- rep('-', length(temp_map[[1]]))
        temp_cancer <- rep('-', length(temp_map[[1]]))

        for(jj in 1:length(tempz1[[1]])){ ## for each splicing event

            tempz2 <- tempz1[jj, ]
            skipped_exons <- tempz2$exons

            for(i in 1:length(skipped_exons)){ ## for each skipped exon, store the ES event ID

                tcga_temp <- tcga_map[tcga_map$Symbol == all_symbol[j], ]
                tcga_temp <- tcga_temp[tcga_temp$Exon == skipped_exons[i], ]
                temp_strand <- tcga_temp$Strand

                for(h in 1:length(nt_start)){ ## for each dbd start
                    ## dbd start and end should be within tcga_temp start and end, in order to consider the dbd to be alternatively spliced

                    if(temp_strand == '-'){
                        wh1 <- (nt_start[h] >= tcga_temp$Chr_Stop)
                        wh2 <- (nt_end[h] <= tcga_temp$Chr_Start)
                    }else{
                        wh1 <- (nt_start[h] >= tcga_temp$Chr_Start)
                        wh2 <- (nt_end[h] <= tcga_temp$Chr_Stop)
                    }
                    
                    if(wh1 & wh2){ 
                        ##-- store the event id and the cancer type ---
                        as_event[h] <- tempz2$as_id
                        # print(paste0(j,':',i))
                    }
                }

            }
        }

        temp_map$AT <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    for(j in 1:length(rest_symbols)){

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == rest_symbols[j])])]

        temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))
        as_event <- rep('-', length(temp_map[[1]]))
        temp_map$AT <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    cat('Cancer',k,'of',length(all_cancer),'done\n')
}

###--- ALTERNATIVE DONOR ---------------------------------------
for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    tempy <- tempx[which(toupper(tempx$symbol) %in% tfs$Gene_Symbol), ]
    tempz <- tempy[tempy$splice_type == 'AD', ]
    all_symbolx <- gtools::mixedsort(unique(tempz$symbol))

    ##--- only the processed TFs --- 1-to-1 mapped transcripts to TFs
    tid <- unlist(lapply(strsplit(unlist(lapply(strsplit(basename(all_filesx), '[_]'),'[[',1)), '[.]'),'[[',1))
    eid <- tf_ensemb_map[which(tf_ensemb_map$Uniprotswissprot %in% tid),]$Ensembl_gene_id
    back_symbols <- tfs[which(tfs$Ensembl_Gene_ID %in% eid), ]$Gene_Symbol 
    ## Out of 1639 curated TFs, 1460 are also present in the exon mapping file
    ### So, the splice event information will be stored for these 1460 TFs
    
    all_symbol <- intersect(all_symbolx, back_symbols)
    rest_symbols <- setdiff(back_symbols, all_symbol)

    for(j in 1:length(all_symbol)){ ## for k=1, you can check for j=6 and i=2

        tempz1 <- tempz[tempz$symbol == all_symbol[j], ]

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == all_symbol[j])])]

        temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))

        nt_start <- temp_map$NT1 ## for each AA position
        nt_end <- temp_map$NT3
        ##------------------------------------------------------------------------
        as_event <- rep('-', length(temp_map[[1]]))
        temp_cancer <- rep('-', length(temp_map[[1]]))

        for(jj in 1:length(tempz1[[1]])){ ## for each splicing event

            tempz2 <- tempz1[jj, ]
            skipped_exons <- tempz2$exons

            for(i in 1:length(skipped_exons)){ ## for each skipped exon, store the ES event ID

                tcga_temp <- tcga_map[tcga_map$Symbol == all_symbol[j], ]
                tcga_temp <- tcga_temp[tcga_temp$Exon == skipped_exons[i], ]
                temp_strand <- tcga_temp$Strand

                for(h in 1:length(nt_start)){ ## for each dbd start
                    ## dbd start and end should be within tcga_temp start and end, in order to consider the dbd to be alternatively spliced

                    if(temp_strand == '-'){
                        wh1 <- (nt_start[h] >= tcga_temp$Chr_Stop)
                        wh2 <- (nt_end[h] <= tcga_temp$Chr_Start)
                    }else{
                        wh1 <- (nt_start[h] >= tcga_temp$Chr_Start)
                        wh2 <- (nt_end[h] <= tcga_temp$Chr_Stop)
                    }
                    
                    if(wh1 & wh2){ 
                        ##-- store the event id and the cancer type ---
                        as_event[h] <- tempz2$as_id
                        # print(paste0(j,':',i))
                    }
                }

            }
        }

        temp_map$AD <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    for(j in 1:length(rest_symbols)){

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == rest_symbols[j])])]

        temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))
        as_event <- rep('-', length(temp_map[[1]]))
        temp_map$AD <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    cat('Cancer',k,'of',length(all_cancer),'done\n')
}


###--- ALTERNATIVE ACCEPTOR ---------------------------------------
for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    tempy <- tempx[which(toupper(tempx$symbol) %in% tfs$Gene_Symbol), ]
    tempz <- tempy[tempy$splice_type == 'AA', ]
    all_symbolx <- gtools::mixedsort(unique(tempz$symbol))

    ##--- only the processed TFs --- 1-to-1 mapped transcripts to TFs
    tid <- unlist(lapply(strsplit(unlist(lapply(strsplit(basename(all_filesx), '[_]'),'[[',1)), '[.]'),'[[',1))
    eid <- tf_ensemb_map[which(tf_ensemb_map$Uniprotswissprot %in% tid),]$Ensembl_gene_id
    back_symbols <- tfs[which(tfs$Ensembl_Gene_ID %in% eid), ]$Gene_Symbol 
    ## Out of 1639 curated TFs, 1460 are also present in the exon mapping file
    ### So, the splice event information will be stored for these 1460 TFs
    
    all_symbol <- intersect(all_symbolx, back_symbols)
    rest_symbols <- setdiff(back_symbols, all_symbol)

    for(j in 1:length(all_symbol)){ ## for k=1, you can check for j=6 and i=2

        tempz1 <- tempz[tempz$symbol == all_symbol[j], ]

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == all_symbol[j])])]

        temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))

        nt_start <- temp_map$NT1 ## for each AA position
        nt_end <- temp_map$NT3
        ##------------------------------------------------------------------------
        as_event <- rep('-', length(temp_map[[1]]))
        temp_cancer <- rep('-', length(temp_map[[1]]))

        for(jj in 1:length(tempz1[[1]])){ ## for each splicing event

            tempz2 <- tempz1[jj, ]
            skipped_exons <- tempz2$exons

            for(i in 1:length(skipped_exons)){ ## for each skipped exon, store the ES event ID

                tcga_temp <- tcga_map[tcga_map$Symbol == all_symbol[j], ]
                tcga_temp <- tcga_temp[tcga_temp$Exon == skipped_exons[i], ]
                temp_strand <- tcga_temp$Strand

                for(h in 1:length(nt_start)){ ## for each dbd start
                    ## dbd start and end should be within tcga_temp start and end, in order to consider the dbd to be alternatively spliced

                    if(temp_strand == '-'){
                        wh1 <- (nt_start[h] >= tcga_temp$Chr_Stop)
                        wh2 <- (nt_end[h] <= tcga_temp$Chr_Start)
                    }else{
                        wh1 <- (nt_start[h] >= tcga_temp$Chr_Start)
                        wh2 <- (nt_end[h] <= tcga_temp$Chr_Stop)
                    }
                    
                    if(wh1 & wh2){ 
                        ##-- store the event id and the cancer type ---
                        as_event[h] <- tempz2$as_id
                        # print(paste0(j,':',i))
                    }
                }

            }
        }

        temp_map$AA <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    for(j in 1:length(rest_symbols)){

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == rest_symbols[j])])]

        temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))
        as_event <- rep('-', length(temp_map[[1]]))
        temp_map$AA <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    cat('Cancer',k,'of',length(all_cancer),'done\n')
}


###--- MUTUALLY EXCLUSIVE ---------------------------------------
for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    tempy <- tempx[which(toupper(tempx$symbol) %in% tfs$Gene_Symbol), ]
    tempz <- tempy[tempy$splice_type == 'AA', ]
    all_symbolx <- gtools::mixedsort(unique(tempz$symbol))

    ##--- only the processed TFs --- 1-to-1 mapped transcripts to TFs
    tid <- unlist(lapply(strsplit(unlist(lapply(strsplit(basename(all_filesx), '[_]'),'[[',1)), '[.]'),'[[',1))
    eid <- tf_ensemb_map[which(tf_ensemb_map$Uniprotswissprot %in% tid),]$Ensembl_gene_id
    back_symbols <- tfs[which(tfs$Ensembl_Gene_ID %in% eid), ]$Gene_Symbol 
    ## Out of 1639 curated TFs, 1460 are also present in the exon mapping file
    ### So, the splice event information will be stored for these 1460 TFs
    
    all_symbol <- intersect(all_symbolx, back_symbols)
    rest_symbols <- setdiff(back_symbols, all_symbol)

    for(j in 1:length(all_symbol)){ ## for k=1, you can check for j=6 and i=2

        tempz1 <- tempz[tempz$symbol == all_symbol[j], ]

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == all_symbol[j])])]

        temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))

        nt_start <- temp_map$NT1 ## for each AA position
        nt_end <- temp_map$NT3
        ##------------------------------------------------------------------------
        as_event <- rep('-', length(temp_map[[1]]))
        temp_cancer <- rep('-', length(temp_map[[1]]))

        for(jj in 1:length(tempz1[[1]])){ ## for each splicing event

            tempz2 <- tempz1[jj, ]
            skipped_exons <- tempz2$exons

            for(i in 1:length(skipped_exons)){ ## for each skipped exon, store the ES event ID

                tcga_temp <- tcga_map[tcga_map$Symbol == all_symbol[j], ]
                tcga_temp <- tcga_temp[tcga_temp$Exon == skipped_exons[i], ]
                temp_strand <- tcga_temp$Strand

                for(h in 1:length(nt_start)){ ## for each dbd start
                    ## dbd start and end should be within tcga_temp start and end, in order to consider the dbd to be alternatively spliced

                    if(temp_strand == '-'){
                        wh1 <- (nt_start[h] >= tcga_temp$Chr_Stop)
                        wh2 <- (nt_end[h] <= tcga_temp$Chr_Start)
                    }else{
                        wh1 <- (nt_start[h] >= tcga_temp$Chr_Start)
                        wh2 <- (nt_end[h] <= tcga_temp$Chr_Stop)
                    }
                    
                    if(wh1 & wh2){ 
                        ##-- store the event id and the cancer type ---
                        as_event[h] <- tempz2$as_id
                        # print(paste0(j,':',i))
                    }
                }

            }
        }

        temp_map$AA <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    for(j in 1:length(rest_symbols)){

        ##-- read the uniprot ensembl map file of the gene --------
        temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
        which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
        which(tfs$Gene_Symbol == rest_symbols[j])])]

        temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))
        as_event <- rep('-', length(temp_map[[1]]))
        temp_map$AA <- as_event
        data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    }

    cat('Cancer',k,'of',length(all_cancer),'done\n')
}











##--- Number of TFs with DBDs affected by ES/AP splicing events -------------------------------
tfs_curated <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
## 1340 out of the 1639 curated TFs have at least one known DBD information ---

##--- store TFs wth DBD info from uniprot ---
TFs_with_DBD <- c()
all_filesx <- list.files(input_dirx, pattern='*.txt', full.names=TRUE)

for(k in 1:length(all_filesx)){

    temp <- data.table::fread(all_filesx[k])
    wh <- which(temp$DBD != '-')
    if(length(wh) != 0){

        tid <- strsplit(basename(all_filesx[k]), '[.]')[[1]][1]
        eid <- tf_ensemb_map[which(tf_ensemb_map$Uniprotswissprot == tid),]$Ensembl_gene_id
        TFs_with_DBD <- union(TFs_with_DBD, tfs_curated[which(tfs$Ensembl_Gene_ID == eid), ]$Gene_Symbol)
    }

}


##------------------------------------------
num_tfs <- c()
num_tf_events <- c()
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI.txt', full.names=TRUE))

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]
    tempy <- tempx[which(toupper(tempx$symbol) %in% tfs$Gene_Symbol), ]
    tempz <- tempy[tempy$splice_type %in% c('ES','AP'), ]
    all_tfs <- unique(tempz$symbol)

    ##-- which of the Tfs have at least one DBD info in the uniprot ----
    cancer_tfs <- intersect(all_tfs, TFs_with_DBD)

    num_tfs <- c(num_tfs, length(cancer_tfs))
    num_tf_events <- c(num_tf_events, length(tempz$as_id))
    
}


##--- plot the number of TFs -----------------
pdata <- data.frame(cancer=all_cancer, count=num_tfs)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of TFs affected by a \n significant splicing event with DBDs", limits=c(0,150)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Affected_TFs_with_DBDs.png"),width=3.5, height=3, dpi=400)

##--- plot the number of splicing events affecting TFs -----------------
pdata <- data.frame(cancer=all_cancer, count=num_tf_events)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of significant events affecting a \n TF with DBDs", limits=c(0,300)) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/Events_affecting_TFs_with_DBDs.png"),width=3.5, height=3, dpi=400)


##----------- Plot the number of TFs where the ES event removes a DBD -------

num_TFs_with_affected_DBD_pt <- c()
num_TFs_with_affected_DBD_ES <- c()
num_TFs_with_affected_DBD_AP <- c()
num_events_affecting_TFs_with_DBD <- c()
TFs_with_affected_DBD <- list()

for(k in 1:length(all_cancer)){

    all_files <- list.files(output_dir, pattern=paste0('*',all_cancer[k],'.txt'), full.names=TRUE)
    count <- 0
    countp <- 0
    all_filesx <- c()
    events <- c()
    for(j in 1:length(all_files)){
        temp <- data.table::fread(all_files[j])
        temp1 <- temp[temp$ES != '-', ]
        temp2 <- temp[temp$AP != '-', ]
        if(nrow(temp1) > 0){
            if(length(which(temp1$DBD != '-')) > 0){
                count <- count+1
                all_filesx <- c(all_filesx, all_files[j])
                events <- c(events, setdiff(unique(temp$ES), '-'))
            }
        }

        if(nrow(temp2) > 0){
            if(length(which(temp2$DBD != '-')) > 0){
                countp <- countp+1
                all_filesx <- c(all_filesx, all_files[j])
                events <- c(events, setdiff(unique(temp$AP), '-'))
            }
        }
    }

    num_TFs_with_affected_DBD_pt <- c(num_TFs_with_affected_DBD_pt, (count/num_tfs[k])*100)
    num_TFs_with_affected_DBD_ES <- c(num_TFs_with_affected_DBD_ES, count)
    num_TFs_with_affected_DBD_AP <- c(num_TFs_with_affected_DBD_AP, countp)
    TFs_with_affected_DBD[[k]] <- unlist(lapply(strsplit(basename(all_filesx), '[_]'),'[[',1))
    num_events_affecting_TFs_with_DBD <- c(num_events_affecting_TFs_with_DBD, length(events))

}


##--- plot the number of splicing events -----------------
pdata <- data.frame(cancer=all_cancer, count=num_TFs_with_affected_DBD_pt, num=num_TFs_with_affected_DBD)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="% of TFs with DBDs affected\n by a significant splicing event", limits=c(0,30)) +
geom_text(aes(label=num), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/TFs_with_affected_DBDs.png"),width=3.5, height=3, dpi=400)


##--- TFs with DBDs affected in multiple cancer types --------------------------------
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
        temp_ovl <- TFs_with_affected_DBD[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- intersect(temp_ovl, TFs_with_affected_DBD[[which(all_cancer == temp_comb[[i]][j])]])
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
scale_y_continuous(name="# of TFs with affected DBD regions", limits=c(0,(max(pdata$count)+5))) +
geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
ggsave(p,filename=paste0(save_dir,"/TFs_DBD_affected_overlap.png"),width=7, height=3, dpi=400)




##----------- Plot the number of TFs where the ES event removes a DBD or an exon participating in forming homodimer -------

num_TFs_with_affected_DBD_pt <- c()
num_TFs_with_affected_DBD <- c()
num_events_affecting_TFs_with_DBD <- c()
TFs_with_affected_DBD <- list()

for(k in 1:length(all_cancer)){

    all_files <- list.files(output_dir, pattern=paste0('*',all_cancer[k],'.txt'), full.names=TRUE)
    count <- 0
    all_filesx <- c()
    events <- c()
    for(j in 1:length(all_files)){
        temp <- data.table::fread(all_files[j])
        temp1 <- temp[temp$ES != '-', ]
        if(nrow(temp1) > 0){
            if(length(which(temp1$DBD != '-')) > 0 | length(which(temp1$HOMODIMER != 0)) > 0){
                count <- count+1
                all_filesx <- c(all_filesx, all_files[j])
                events <- c(events, setdiff(unique(temp$ES), '-'))
            }
        }
    }

    num_TFs_with_affected_DBD_pt <- c(num_TFs_with_affected_DBD_pt, (count/num_tfs[k])*100)
    num_TFs_with_affected_DBD <- c(num_TFs_with_affected_DBD, count)
    TFs_with_affected_DBD[[k]] <- unlist(lapply(strsplit(basename(all_filesx), '[_]'),'[[',1))
    num_events_affecting_TFs_with_DBD <- c(num_events_affecting_TFs_with_DBD, length(events))

}


##--- plot the number of splicing events -----------------
pdata <- data.frame(cancer=all_cancer, count=num_TFs_with_affected_DBD_pt, num=num_TFs_with_affected_DBD)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity",position=position_dodge())+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="% of TFs with DBDs affected\n by a significant splicing event", limits=c(0,30)) +
geom_text(aes(label=num), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
# scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/TFs_with_affected_DBDs_dimer.png"),width=3.5, height=3, dpi=400)







