##############################################################################################
# Purpose: create map of TFs and exons and DBDs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)

tfs <- data.table::fread('../data/filtered_TFs.txt', sep='\t')

## Proteins from Uniprot...these are reviewed entries --- downloaded 5th March 2025 -------------------
system("wget -O ../../../public_data/uniprot_sprot.dat.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz")
system("gunzip --force ../../../public_data/uniprot_sprot.dat.gz")

## save separate files for uniprot entries ---- takes time and hence should be run once ---
uniprot_dir <- '../../../public_data/UniProt_database'
unlink(uniprot_dir,recursive=TRUE)
dir.create(uniprot_dir)
system("awk '/^ID/{if(NR!=1){for(i=0;i<j;i++)print a[i]>\"../../../public_data/UniProt_database/file\"k;j=0;k++;}a[j++]=$0;next}{a[j++]=$0;}END{for(i=0;i<j;i++)print a[i]>\"../../../public_data/UniProt_database\"k}' i=0 k=1  ../../../public_data/uniprot_sprot.dat")


## filter to only retain Homo sapiens -----
allfiles <- list.files(uniprot_dir, full.names=TRUE)
output_dir <- '../../../public_data/UniProt_database_HS'
## unlink(output_dir, recursive=TRUE)
dir.create(output_dir)
counter <- 0
for(k in 1:length(allfiles)){

    tempf <- readLines(allfiles[k])

    for(j in 1:length(tempf)){

        temp <- tempf[j]
        temp2 <- substr(temp, 1,2)

        if((temp2 == 'OS') & (temp %like% 'Homo sapiens')){
            tf <- strsplit(tempf[1],'\\s+')
            filename <- paste0(output_dir,'/',tf[[1]][2],'_',strsplit(tf[[1]][3],';')[[1]][1],'.txt')
            file.copy(allfiles[k], filename)
            counter <- counter+1                            
            break
        }

    }

    cat(k,' of ', length(allfiles), ' done\n')

}


##----- save uniprot ensembl mapping -------------------------------
output_dir <- '../../../public_data/UniProt_database_HS'
allfiles <- list.files(output_dir, full.names=TRUE)
temp_trans <- c()
temp_prot <- c()
temp_gene <- c()
temp_unip <- c()
temp_sym <- c()
temp_hgnc <- c()

for(k in 1:length(allfiles)){
    tempf <- readLines(allfiles[k])
    temp2 <- substr(tempf, 1,2)
    whp <- which(temp2 == 'DR')
    tempf1 <- tempf[whp]
    temp3 <- substr(tempf1, 6,12)
    whp <- which(temp3 == 'Ensembl')
    if(length(whp) != 0){
        tempf2 <- tempf1[whp]
        tempf3 <- strsplit(tempf2, '[;]')
        temp_trans <- c(temp_trans, trimws(unlist(lapply(tempf3,'[[',2))) )
        temp_prot <- c(temp_prot, trimws(unlist(lapply(tempf3,'[[',3))) )
        tempf4 <- strsplit(unlist(lapply(tempf3,'[[',4)), '[\\[]')
        if((length(tempf4) > 1) & min(lengths(tempf4)) > 1){
            for(j in 1:length(tempf4)){
                ttx <- trimws(tempf4[[j]][1])
                ttu <- trimws(tempf4[[j]][2])
                temp_gene <- c(temp_gene, substr(ttx,1,nchar(ttx)-1))
                temp_unip <- c(temp_unip, substr(ttu,1,nchar(ttu)-1))
            }
        }else{
            for(j in 1:length(tempf4)){
                ttx <- trimws(tempf4[[j]][1])
                temp_gene <- c(temp_gene, substr(ttx,1,nchar(ttx)-1))
                wha <- which(temp2 == 'AC')
                # take the canonical uniprot identifier
                temp1 <- unlist(lapply(strsplit(tempf[wha[1]], '[;]'), '[[',1))
                temp11 <- lapply(strsplit(temp1,'\\s+'), '[[',2)
                temp_unix <- paste0(temp11[[1]],'-1')
                temp_unip <- c(temp_unip, temp_unix)
            }
        }

        temp3 <- substr(tempf1, 6,9)
        whp <- which(temp3 == 'HGNC')
        if(length(whp) != 0){
            tty <- trimws(strsplit(tempf1[whp], '[;]')[[1]][3])
            temp_hgnc <- c(temp_hgnc, rep(substr(tty,1,nchar(tty)-1), length(tempf4)))
        }else{
            temp_hgnc <- c(temp_hgnc, rep('', length(tempf4)))
        }
    }
    cat(k,' of ', length(allfiles), ' done\n')
}

map_data <- data.frame(Uniprotswissprot=temp_unip, Ensembl_gene_id=temp_gene, 
    Ensembl_protein_id=temp_prot, Ensembl_transcript_id=temp_trans, HGNC_symbol=temp_hgnc)

data.table::fwrite(map_data, '../data/Uniprot_ensembl_map.txt', sep='\t', row.names=FALSE, quote=FALSE)

##------------------------------------------------------------------

### saving uniprot AA sequences ------------------------------------
## right now saving the canonical sequence of a protein. Other transcripts could be longer
## Next version should consider saving the longest coding seqeunce---
uniprot_HS_dir <- '../../../public_data/UniProt_database_HS'
allfiles <- list.files(uniprot_HS_dir, full.names=TRUE)
output_file <- '../data/uniprot_sequences_HS.txt'
if(file.exists(output_file)){file.remove(output_file)}
uniprot_name <- c()
uniprot_id <- c()

# save all sequences in a vector
uniseq_aa <- list()
for(k in 1:length(allfiles)){

    temp <- readLines(allfiles[k])
    tempseq <- substr(temp, 1,2)

    uniprot_name <- c(uniprot_name, strsplit(basename(allfiles[k]), '[.]')[[1]][1])

    wha <- which(tempseq == 'AC')
    # take the canonical uniprot identifier
    temp1 <- unlist(lapply(strsplit(temp[wha[1]], '[;]'), '[[',1))
    temp11 <- lapply(strsplit(temp1,'\\s+'), '[[',2)
    temp_uniprot <- temp11[[1]]

    uniprot_id <- c(uniprot_id, temp_uniprot)

    wh1 <- which(tempseq == 'SQ')+1
    wh2 <- which(tempseq == '//')-1

    temps1 <- gsub(' ','',paste(temp[wh1:wh2], collapse=''))
    seqinr::write.fasta(temps1, temp_uniprot, output_file, open='a', nbchar=60, as.string=FALSE)
    cat(k,' of ', length(allfiles), ' done\n')

}

uniprot_name_id_map <- data.frame(Name=uniprot_name, ID=uniprot_id)
data.table::fwrite(uniprot_name_id_map,'../data/uniprot_name_id_map.txt',sep='\t', row.names=FALSE, quote=FALSE)

##---Ensembl mapping -------------------------------------------------------
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

cleanSeq <- function(x){
    last_flag <- substrRight(x, 1)
    if(last_flag == "*"){
        slen <- nchar(x)-1
        x <- substr(x, 1, slen)
    }
    return(x)
}

##--- TCGASPLiceseq uses GRCh37 build ------
# Ensemble genome download -------------------------
xx <- "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
system(paste0("wget -O ../../../public_data/",xx," https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/",xx))
system(paste0("gunzip ../../../public_data/",xx))

# download GTF file from Ensembl
xx <- "Homo_sapiens.GRCh38.113.gtf.gz"
system(paste0("wget -O ../../../public_data/",xx," http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/",xx))
system(paste0("gunzip ../../../public_data/",xx))

xx <- "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
system(paste0("wget -O ../../../public_data/",xx," https://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/",xx))
system(paste0("gunzip ../../../public_data/",xx))

# https://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
# download GTF file from Ensembl
xx <- "Homo_sapiens.GRCh37.87.gtf.gz"
system(paste0("wget -O ../../../public_data/",xx," http://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/",xx))
system(paste0("gunzip ../../../public_data/",xx))

# preprocess
system('sh ./preprocess_exons.sh')
exons <- fread('../data/processed_exons.tmp',sep='\t')
cds <- fread('../data/processed_cds.tmp',sep='\t')


xxs <- strsplit(exons$V9, '[;]')
# gene_id
wh1 <- lapply(xxs, function(x) which(x %like% 'gene_id'))
gene_id1 <- mapply(function(x, y) x[y], xxs, wh1)
gene_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(gene_id1),'\\s+'), '[[', 2)))

# transcript_id
wh1 <- lapply(xxs, function(x) which(x %like% 'transcript_id'))
transcript_id1 <- mapply(function(x, y) x[y], xxs, wh1)
transcript_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(transcript_id1),'\\s+'), '[[', 2)))

# exon_number
wh1 <- lapply(xxs, function(x) which(x %like% 'exon_number'))
exon_number1 <- mapply(function(x, y) x[y], xxs, wh1)
exon_number <- gsub('\\"', "", unlist(lapply(strsplit(trimws(exon_number1),'\\s+'), '[[', 2)))

# transcript_biotype
wh1 <- lapply(xxs, function(x) which(x %like% 'transcript_biotype'))
transcript_biotype1 <- mapply(function(x, y) x[y], xxs, wh1)
transcript_biotype <- gsub('\\"', "", unlist(lapply(strsplit(trimws(transcript_biotype1),'\\s+'), '[[', 2)))

# exon_id
wh1 <- lapply(xxs, function(x) which(x %like% 'exon_id'))
exon_id1 <- mapply(function(x, y) x[y], xxs, wh1)
exon_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(exon_id1),'\\s+'), '[[', 2)))

exons$gene_id <- gene_id
exons$transcript_id <- transcript_id
exons$transcript_biotype <- transcript_biotype
exons$exon_number <- exon_number
exons$exon_id <- exon_id

exons_f <- exons[,-9]
data.table::fwrite(exons_f, '../data/Ensembl_exons.txt', sep='\t', quote=FALSE, row.names=FALSE)
exons_f <- data.table::fread('../data/Ensembl_exons.txt', header=TRUE)


#### CDS
xxs <- strsplit(cds$V9, '[;]')
# gene_id
wh1 <- lapply(xxs, function(x) which(x %like% 'gene_id'))
gene_id1 <- mapply(function(x, y) x[y], xxs, wh1)
gene_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(gene_id1),'\\s+'), '[[', 2)))

# transcript_id
wh1 <- lapply(xxs, function(x) which(x %like% 'transcript_id'))
transcript_id1 <- mapply(function(x, y) x[y], xxs, wh1)
transcript_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(transcript_id1),'\\s+'), '[[', 2)))

# exon_number
wh1 <- lapply(xxs, function(x) which(x %like% 'exon_number'))
exon_number1 <- mapply(function(x, y) x[y], xxs, wh1)
exon_number <- gsub('\\"', "", unlist(lapply(strsplit(trimws(exon_number1),'\\s+'), '[[', 2)))

# transcript_biotype
wh1 <- lapply(xxs, function(x) which(x %like% 'transcript_biotype'))
transcript_biotype1 <- mapply(function(x, y) x[y], xxs, wh1)
transcript_biotype <- gsub('\\"', "", unlist(lapply(strsplit(trimws(transcript_biotype1),'\\s+'), '[[', 2)))

cds$gene_id <- gene_id
cds$transcript_id <- transcript_id
cds$transcript_biotype <- transcript_biotype
cds$exon_number <- exon_number
cds$nt_len <- (cds$V5-cds$V4)+1

cds_f <- cds[,-9]
fwrite(cds_f, '../data/Ensembl_exon_cds.txt', sep='\t', quote=FALSE, row.names=FALSE)
cds_f <- fread('../data/Ensembl_exon_cds.txt', header=TRUE)


# # get the uniprot emsemble mapping ##### USE OF BIOMART --------------------------------------------------------
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror='useast') # useast uswest asia

##-- read the genome -----
# priassm <- '../../../public_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
priassm <- '../../../public_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa'
genome <- seqinr::read.fasta(priassm)

all_map <- data.table::fread('../data/Uniprot_ensembl_map.txt')
all_map$Ensembl_gene_id <- unlist(lapply(strsplit(all_map$Ensembl_gene_id, '[.]'),'[[',1))
temp_map <- all_map[all_map$Ensembl_gene_id %in% tfs[[1]], ]
temp_map$Uniprotswissprot <- unlist(lapply(strsplit(temp_map$Uniprotswissprot, '[-]'),'[[',1))
temp_map$Ensembl_transcript_id <- unlist(lapply(strsplit(temp_map$Ensembl_transcript_id, '[.]'),'[[',1))

# temp_map <- temp_map[temp_map$UNIPROT == '', ]
# temp_map1 <- unique(temp_map[,c(1,13)])
# temp_count <- plyr::count(temp_map1$Uniprotswissprot)
# wh1 <- which(temp_map$Uniprotswissprot %in% temp_count[[1]][which(temp_count$freq > 1)])
# temp_count <- plyr::count(temp_map1$Ensembl_gene_id)
# wh2 <- which(temp_map$Ensembl_gene_id %in% temp_count[[1]][which(temp_count$freq > 1)])
# wh <- union(wh1, wh2)
# temp_map <- temp_map[-wh, ]
## 2720 out of the 2765 TFs have one-to-one Ensembl ids and uniprot ids mapped


##--- select the transcript corresponding to the sequence present in the uniprot database ---
##--- using global alignment ----
uniprot_seqs <- seqinr::read.fasta('../data/uniprot_sequences_HS.txt', seqtype="AA", whole.header=TRUE)
#unique(unlist(lapply(strsplit(temp_map$UNIPROT,'[-]'),'[[',1)))
ensemblids <- unique(temp_map$Ensembl_gene_id)


# countl <- rep(0, length(ensemblids))
count_gap <- rep(0, length(ensemblids))
iden <- rep(0,length(ensemblids))
bh_transcripts <- rep('',length(ensemblids))
uniprotids <- c()
for(k in 1:length(ensemblids)){

    alldata <- temp_map[temp_map$Ensembl_gene_id == ensemblids[k], ]

    temp_uni_id <- unique(unlist(lapply(strsplit(alldata$Uniprotswissprot,'[-]'),'[[',1)))[1]
    uniprotids <- c(uniprotids, temp_uni_id)
    uni_seq <- seqinr::getSequence(uniprot_seqs[temp_uni_id])
    uni_seq_len <- seqinr::getLength(uniprot_seqs[temp_uni_id])
    toalign1 <- toupper(paste(uni_seq[[1]],collapse=''))

    all_transcripts <- unlist(lapply(strsplit(unique(alldata$Ensembl_transcript_id),'[.]'),'[[',1))

    if(!identical(all_transcripts, character(0))){

        for(j in 1:length(all_transcripts)){

            temp_cds <- cds_f[cds_f$transcript_id == all_transcripts[j], ]
            
            if(nrow(temp_cds) > 0){
                dna_strand <- temp_cds$V7[[1]][1]
                temp_sum <- sum(temp_cds$nt_len)/3

                if(temp_sum == uni_seq_len){

                    temps <- seqinr::getSequence(genome[temp_cds[[1]][1]])[[1]]

                    allseq <- c()
                    for(i in 1:length(temp_cds[[1]])){
                        start <- temp_cds$V4[i]
                        end <- temp_cds$V5[i]
                        tempseq <- temps[start:end]
                        if(dna_strand == "-"){
                            tempseq <- rev(seqinr::comp(tempseq))
                        }
                        allseq <- c(allseq, tempseq)
                    }
                    
                    toalign2 <- toupper(paste(seqinr::translate(allseq),collapse=''))

                    xxx <- Biostrings::pairwiseAlignment(pattern=toalign1, subject=toalign2,
                                 substitutionMatrix = "BLOSUM62", 
                                 gapOpening = -2,
                                 gapExtension = -8, 
                                 scoreOnly = FALSE)

                    iso1 <- as.vector(stringr::str_split_fixed(as.character(Biostrings::pattern(xxx)), pattern='', n=nchar(as.character(Biostrings::pattern(xxx)))))
                    iso2 <- as.vector(stringr::str_split_fixed(as.character(Biostrings::subject(xxx)), pattern='', n=nchar(as.character(Biostrings::subject(xxx)))))

                    gap1 <- length(which(iso1 == '-')) # number of gaps in first seq
                    gap2 <- length(which(iso2 == '-')) # number of gaps in second seq
                    gap <- gap1+gap2
                    count_gap[k] <- gap
                    iden[k] <- 100-Biostrings::pid(xxx, type='PID3') ##
                    if(iden[k] == 0 && gap == 0){
                        bh_transcripts[k] <- all_transcripts[j]
                        break
                    }
                }
            }
        }
    }
    cat('Protein',k,' of ', length(ensemblids), ' done\n')
}

temp_mapx <- data.frame(Uniprotswissprot=uniprotids, Ensembl_transcript_id=bh_transcripts)
temp_mapx <- temp_mapx[temp_mapx$Ensembl_transcript_id != '',]
temp_mapz <- merge(temp_mapx, temp_map, by='Ensembl_transcript_id')
# temp_mapz <- unique(temp_mapy[, c(1,2,3,14)])
# wh <- which(temp_mapz$Ensembl_transcript == '')
# temp_mapz <- temp_mapz[-wh,]
## out of 2720 TFs, 2474 have a 100% transcipt match between uniprot protein and an Ensembl transcript
## intersect between uniprots here and TFs -- intersect(temp_mapz$Ensembl_gene_id, tfs[[1]])
temp_mapz <- temp_mapz[,-2]
colnames(temp_mapz) <- c('Ensembl_transcript_id','Uniprotswissprot',colnames(temp_mapz)[-c(1,2)])
data.table::fwrite(temp_mapz,'../data/TF_ensembl_uniprot.txt', sep='\t', row.names=FALSE, quote=FALSE)


##------- do the exon mapping to uniprot sequences -------------------------------------
temp_map <- data.table::fread('../data/TF_ensembl_uniprot.txt')
uniprot_uni_ids <- temp_map$Uniprotswissprot
transcript_uni_ids <- temp_map$Ensembl_transcript_id

##-- map exons to uniprot sequences -------------
##-- store the mappings -------------------------
store_dir <- '../data/uniprot_Ensembl_Exon_map'
if(dir.exists(store_dir)){
    unlink(store_dir, recursive=TRUE)
}
dir.create(store_dir)

for(k in 1:length(transcript_uni_ids)){

    # get the uniprot sequence
    uniprot_seq <- seqinr::getSequence(uniprot_seqs[uniprot_uni_ids[k]])[[1]]

    # get all exons for this transcript
    allexons <- exons_f[exons_f$transcript_id == transcript_uni_ids[k], ]

    # get exon numbers from cds...denoting the coding exons for this transcript
    cexon_num <- cds_f[cds_f$transcript_id == transcript_uni_ids[k], ]

    if(nrow(cexon_num) != 0){ ## if the concerned transcript is not in the CDS file from Ensembl

        cexon_num <- cexon_num[order(cexon_num$exon_number), ]

        # get exons
        cexon <- allexons[allexons$exon_number %in% cexon_num$exon_number, ]

        # sort exons by genomic start
        cexon <- cexon[order(cexon$exon_number), ]

        # for each exon
        exon_entry <- rep('',length(uniprot_seq))
        exon_num_entry <- rep('',length(uniprot_seq))
        chrmsm <- rep('',length(uniprot_seq))
        strand_dr <- rep('',length(uniprot_seq))
        startp <- rep('',length(uniprot_seq))
        endp <- rep('',length(uniprot_seq))
        startpp <- rep('',length(uniprot_seq))
        endpp <- rep('',length(uniprot_seq))
        start_pos <- 1
        end_pos <- 0
        previous <- 0
        genome_pos <- c()

        for(j in 1:length(cexon[[1]])){

            ## genomic positions
            genome_pos <- c(genome_pos, seq(cexon_num$V4[j],cexon_num$V5[j]))

            temp_pos <- (cexon_num$nt_len[j]+previous)/3

            if(temp_pos != 0){

                temp_pos1 <- floor(temp_pos)
                diff <- temp_pos-temp_pos1

                if(diff == 0){ # integer temp_pos
                    end_pos <- end_pos+temp_pos1
                    exon_entry[start_pos:end_pos] <- cexon$exon_id[j]
                    exon_num_entry[start_pos:end_pos] <- cexon$exon_number[j]
                    chrmsm[start_pos:end_pos] <- cexon$V1[j]
                    strand_dr[start_pos:end_pos] <- cexon$V7[j]
                    startp[start_pos:end_pos] <- cexon$V4[j]
                    endp[start_pos:end_pos] <- cexon$V5[j]
                    startpp[start_pos:end_pos] <- cexon_num$V4[j]
                    endpp[start_pos:end_pos] <- cexon_num$V5[j]
                    start_pos <- end_pos+1
                    previous <- 0
                }else if(diff < 0.5){
                    end_pos <- end_pos+temp_pos1
                    exon_entry[start_pos:end_pos] <- cexon$exon_id[j]
                    exon_num_entry[start_pos:end_pos] <- cexon$exon_number[j]
                    chrmsm[start_pos:end_pos] <- cexon$V1[j]
                    strand_dr[start_pos:end_pos] <- cexon$V7[j]
                    startp[start_pos:end_pos] <- cexon$V4[j]
                    endp[start_pos:end_pos] <- cexon$V5[j]
                    startpp[start_pos:end_pos] <- cexon_num$V4[j]
                    endpp[start_pos:end_pos] <- cexon_num$V5[j]
                    start_pos <- end_pos+1
                    previous <- 1
                }else{
                    end_pos <- end_pos+temp_pos1+1
                    exon_entry[start_pos:end_pos] <- cexon$exon_id[j]
                    exon_num_entry[start_pos:end_pos] <- cexon$exon_number[j]
                    chrmsm[start_pos:end_pos] <- cexon$V1[j]
                    strand_dr[start_pos:end_pos] <- cexon$V7[j]
                    startp[start_pos:end_pos] <- cexon$V4[j]
                    endp[start_pos:end_pos] <- cexon$V5[j]
                    startpp[start_pos:end_pos] <- cexon_num$V4[j]
                    endpp[start_pos:end_pos] <- cexon_num$V5[j]
                    start_pos <- end_pos+1
                    previous <- -1
                }

            }
        
        }

        nt1 <- c()
        nt2 <- c()
        nt3 <- c()
        for(i in seq(from=1, to=length(genome_pos), by=3)){
            nt1 <- c(nt1, genome_pos[i])
            nt2 <- c(nt2, genome_pos[i+1])
            nt3 <- c(nt3, genome_pos[i+2])
        }

        Data <- data.frame(UNIPROT_SEQ_NUM=seq(1,length(uniprot_seq)), UNIPROT=uniprot_seq, EXON=exon_entry, 
            EXON_NUM=exon_num_entry, CHR_NUM=chrmsm, CHR_START_POS=startp, CHR_END_POS=endp, 
            CHR_START_POS_CODING=startpp, CHR_END_POS_CODING=endpp, CHR_STAND=strand_dr, 
            NT1=nt1,NT2=nt2,NT3=nt3)
        
        fwrite(Data, paste0(store_dir,'/',uniprot_uni_ids[k],'.txt'), row.names=FALSE, sep='\t', quote=FALSE)

    }
    
    cat('Protein',k,' of ', length(transcript_uni_ids), ' done\n')

}



##----------- Add DBD information for the TFs from UniProt files ------------------------
# Check that it doesn't match any non-number
numbers_only <- function(x) !grepl("\\D", x)
store_dir <- '../data/uniprot_Ensembl_Exon_map'
# data.table::fwrite(uniprot_name_id_map,'../data/uniprot_name_id_map.txt',sep='\t', row.names=FALSE, quote=FALSE)

uniprot_name_id_map <- data.table::fread('../data/uniprot_name_id_map.txt')
allfiles <- list.files(store_dir,full.names=TRUE)
DBDs <- c()
pro_dbd <- c()
DBDs_list <- vector(mode = "list", length = length(allfiles))

tprot <- c()
tdbd <- c()
tgene <- c()
tevi <- c()

for(k in 1:length(allfiles)){

    temp_uniprot <- strsplit(basename(allfiles[k]), '[.]')[[1]][1]
    temp_file <- data.table::fread(allfiles[k])

    ##---- read uniprot file ----
    uni_name <- uniprot_name_id_map$Name[which(uniprot_name_id_map$ID == temp_uniprot)] 
    uniprot_file <- readLines(list.files(uniprot_HS_dir, pattern=paste0('^', uni_name), full.names=TRUE))

    temp <- substr(uniprot_file, 1,2)
    temp_uniprot_file <- uniprot_file[which(temp == 'FT')]
    dna_bind_pos1 <- which(temp_uniprot_file %like% 'DNA_BIND')
    dna_bind_pos2 <- which(temp_uniprot_file %like% 'ZN_FING')
    dna_bind_pos <- union(dna_bind_pos1, dna_bind_pos2)
    dna_bind <- rep('-', length(temp_file[[1]]))
    PROSITE <- rep('-', length(temp_file[[1]]))
    EVID <- rep('-', length(temp_file[[1]]))
    temp_binds <- c()
    ## loop for DNA_BIND and ZN_FING---
    if(length(dna_bind_pos) != 0){
        for(j in 1:length(dna_bind_pos)){
            temp_dna <- strsplit(temp_uniprot_file[dna_bind_pos[j]],"\\s+")[[1]][3]
            dna_bseq <- seq(as.numeric(strsplit(temp_dna, '[..]')[[1]][1]), as.numeric(strsplit(temp_dna, '[..]')[[1]][3]))

            temp_next <- temp_uniprot_file[dna_bind_pos[j]+1]
            if(temp_next %like% 'note='){
                dna_bind_domain <- gsub('\\"',"",strsplit(temp_next,"=")[[1]][2])
            }else if(temp_next %like% 'evidence='){
                dna_bind_domain <- 'Unknown'
                evidence <- strsplit(gsub('\\"',"",strsplit(temp_next,"=")[[1]][2]),"[|]")[[1]][1]
                prorule <- strsplit(gsub('\\"',"",strsplit(temp_next,"[|]")[[1]][2]),'[:]')[[1]][2]
            }else{
                dna_bind_domain <- 'Unknown'
                evidence <- NA
                prorule <- NA
            }

            temp_next <- temp_uniprot_file[dna_bind_pos[j]+2]
            if(temp_next %like% 'evidence='){
                evidence <- strsplit(gsub('\\"',"",strsplit(temp_next,"=")[[1]][2]),"[|]")[[1]][1]
                prorule <- strsplit(gsub('\\"',"",strsplit(temp_next,"[|]")[[1]][2]),'[:]')[[1]][2]
            }else{
                evidence <- NA
                prorule <- NA
            }
            
            dna_bind_domain <- strsplit(dna_bind_domain,'[;]')[[1]][1]

            dnaxx <- strsplit(dna_bind_domain,'[ ]')
            dnaxxf <- dnaxx[[1]][length(dnaxx[[1]])]
            loop <- length(dnaxx[[1]])-1
            if(numbers_only(dnaxxf)){
                dnaxxp <- dnaxx[[1]][1]
                if(loop > 1){
                    for(hh in 2:loop){
                        dnaxxp <- paste(dnaxxp,dnaxx[[1]][hh])
                    }
                } 
                dna_bind_domain <- dnaxxp
            }
            
            wh <- which(temp_file$UNIPROT_SEQ_NUM %in% dna_bseq)
            dna_bind[wh] <- dna_bind_domain
            PROSITE[wh] <- substr(prorule,1,8)
            EVID[wh] <- evidence

            DBDs <- c(DBDs, dna_bind_domain)
            pro_dbd <- c(pro_dbd, substr(prorule,1,8))
            tgene <- c(tgene, uni_name)
            tevi <- c(tevi, evidence)
            temp_binds <- c(temp_binds, dna_bind_domain)
        }
    }
    # if('Homeobox' %in% temp_binds){
    #     break
    # }

    ##--- loop for bHLH ------
    to_consider_domains <- c('bHLH','bZIP','MADS-box')
    dna_bind_pos <- which(temp_uniprot_file %like% 'DOMAIN')

    if(length(dna_bind_pos) != 0) {
        for(j in 1:length(dna_bind_pos)){
            # dna_bind_domain <- gsub('\\"',"",strsplit(strsplit(temp_uniprot_file[dna_bind_pos[j]+1],"\\s+")[[1]][2],'=')[[1]][2])
            dna_bind_domain <- gsub('\\"',"",strsplit(temp_uniprot_file[dna_bind_pos[j]+1],"=")[[1]][2])
            dna_bind_domain <- strsplit(dna_bind_domain,'[;]')[[1]][1]
            if(dna_bind_domain %in% to_consider_domains){
                temp_dna <- strsplit(temp_uniprot_file[dna_bind_pos[j]],"\\s+")[[1]][3]
                dna_bseq <- seq(as.numeric(strsplit(temp_dna, '[..]')[[1]][1]), as.numeric(strsplit(temp_dna, '[..]')[[1]][3]))

                temp_next <- temp_uniprot_file[dna_bind_pos[j]+2]
                if(temp_next %like% 'evidence='){
                    evidence <- strsplit(gsub('\\"',"",strsplit(temp_next,"=")[[1]][2]),"[|]")[[1]][1]
                    prorule <- strsplit(gsub('\\"',"",strsplit(temp_next,"[|]")[[1]][2]),'[:]')[[1]][2]
                }else{
                    evidence <- NA
                    prorule <- NA
                }

                wh <- which(temp_file$UNIPROT_SEQ_NUM %in% dna_bseq)
                dna_bind[wh] <- dna_bind_domain
                PROSITE[wh] <- substr(prorule,1,8)
                EVID[wh] <- evidence
            
                DBDs <- c(DBDs, dna_bind_domain)
                pro_dbd <- c(pro_dbd, substr(prorule,1,8))
                tgene <- c(tgene, uni_name)
                tevi <- c(tevi, evidence)
                temp_binds <- c(temp_binds, dna_bind_domain)
            } 
        }
    }
    
    temp_file$DBD <- dna_bind
    # temp_file$PROSITE <- PROSITE
    # temp_file$EVIDENCE <- EVID

    DBDs_list[[k]] <- temp_binds
    data.table::fwrite(temp_file, paste0(store_dir,'/',temp_uniprot,'.txt'), row.names=FALSE, sep='\t', quote=FALSE)

    tprot <- c(tprot, temp_uniprot)
    tdbd <- c(tdbd, paste(unique(setdiff(temp_binds,'-')), collapse=';'))
    cat('TF',k,'of',length(allfiles),'done\n')

}
