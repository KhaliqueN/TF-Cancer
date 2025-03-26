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
input_dirx <- '../data/uniprot_Ensembl_Exon_map'
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
                        whx <- (nt_end[h] <= tcga_temp$Chr_Stop)
                    }else{
                        whx <- (nt_end[h] >= tcga_temp$Chr_Stop)
                    }
                    
                    if(whx){ 
                        ##-- store the event id and the cancer type ---
                        as_event[h] <- tempz2$as_id
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
                        as_event[h] <- tempz2$as_id
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
                        as_event[h] <- tempz2$as_id
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
    tempz <- tempy[tempy$splice_type == 'ME', ]
    
    if(nrow(tempz) != 0){
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
                skipped_exons <- strsplit(tempz2$exons,'[|]')[[1]]

                for(i in 1:length(skipped_exons)){ ## for each skipped exon, store the ES event ID
                    temp_exon <- strsplit(skipped_exons[i],'[:]')[[1]]

                    for(ii in 1:length(temp_exon)){

                        tcga_temp <- tcga_map[tcga_map$Symbol == all_symbol[j], ]
                        tcga_temp <- tcga_temp[tcga_temp$Exon == temp_exon[ii], ]
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
            }

            temp_map$ME <- as_event
            data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

        }

        for(j in 1:length(rest_symbols)){

            ##-- read the uniprot ensembl map file of the gene --------
            temp_uniprot <- tf_ensemb_map$Uniprotswissprot[ 
            which(tf_ensemb_map$Ensembl_gene_id == tfs$Ensembl_Gene_ID[
            which(tfs$Gene_Symbol == rest_symbols[j])])]

            temp_map <- data.table::fread(paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'))
            as_event <- rep('-', length(temp_map[[1]]))
            temp_map$ME <- as_event
            data.table::fwrite(temp_map,paste0(output_dir,'/',temp_uniprot,'_',all_cancer[k],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

        }
    }
    

    cat('Cancer',k,'of',length(all_cancer),'done\n')
}


