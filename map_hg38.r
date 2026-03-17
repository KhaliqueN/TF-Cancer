##############################################################################################
# Purpose: cmap hg38
##############################################################################################

rm(list=ls())

library(data.table)
# library(ggplot2)
# library(GenomicDataCommons)
# library(biomaRt)
# library(seqinr)


# Ensemble genome download -------------------------
xx <- "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
system(paste0("wget -O data/",xx," https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/",xx))
system(paste0("gunzip data/",xx))

# download GTF file from Ensembl
xx <- "Homo_sapiens.GRCh38.113.gtf.gz"
system(paste0("wget -O data/",xx," http://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/",xx))
system(paste0("gunzip data/",xx))


# xx <- "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
# system(paste0("wget -O data/",xx," https://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/",xx))
# system(paste0("gunzip data/",xx))

# # https://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
# # download GTF file from Ensembl
# xx <- "Homo_sapiens.GRCh37.87.gtf.gz"
# system(paste0("wget -O data/",xx," http://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/",xx))
# system(paste0("gunzip data/",xx))

# preprocess
system('sh ./preprocess_exon_hg38.sh')
exons <- fread('data/processed_exons.tmp',sep='\t')
cds <- fread('data/processed_cds.tmp',sep='\t')


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
data.table::fwrite(exons_f, 'data/Ensembl_exons_hg38.txt', sep='\t', quote=FALSE, row.names=FALSE)
exons_f <- data.table::fread('data/Ensembl_exons_hg38.txt', header=TRUE)


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
fwrite(cds_f, 'data/Ensembl_exon_cds_hg38.txt', sep='\t', quote=FALSE, row.names=FALSE)
cds_f <- fread('data/Ensembl_exon_cds_hg38.txt', header=TRUE)



##------- do the exon mapping to uniprot sequences -------------------------------------
temp_map <- data.table::fread('data/TF_ensembl_uniprot.txt')
uniprot_uni_ids <- temp_map$Uniprotswissprot
transcript_uni_ids <- temp_map$Ensembl_transcript_id

##-- are the transcript ids from grch37 present in hg38?
dd <- intersect(transcript_uni_ids, exons_f$transcript_id) ## YES, every transcript id prevails

##-- map exons to uniprot sequences -------------
##-- store the mappings -------------------------
in_dir <- 'data/uniprot_Ensembl_Exon_map_DBD_ED_AS_Lambourne'
all_files <- list.files(in_dir, full.names=TRUE)
store_dir <- 'data/uniprot_Ensembl_Exon_map_DBD_ED_AS_Lambourne_hg38'
if(dir.exists(store_dir)){
    unlink(store_dir, recursive=TRUE)
}
dir.create(store_dir)
processed <- 0
not_processed <- c()

for(k in 1:length(all_files)){

    tempf <- data.table::fread(all_files[k])

    uniprot_seq <- tempf$UNIPROT
    temp_transcript <- temp_map[temp_map$Uniprotswissprot == strsplit(basename(all_files[k]),'[_]')[[1]][1], ]$Ensembl_transcript_id

    # get all exons for this transcript
    allexons <- exons_f[exons_f$transcript_id == temp_transcript, ]

    # get exon numbers from cds...denoting the coding exons for this transcript
    cexon_num <- cds_f[cds_f$transcript_id == temp_transcript, ]

    if(nrow(cexon_num) != 0){ ## if the concerned transcript is present in the CDS file from Ensembl

        processed <- processed+1

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

        if(length(uniprot_seq) != length(nt1)){
            not_processed <- c(not_processed, k)
            next
        }

        Data <- data.frame(UNIPROT_SEQ_NUM_hg38=seq(1,length(uniprot_seq)), UNIPROT_hg38=uniprot_seq, EXON_hg38=exon_entry, 
            EXON_NUM_hg38=exon_num_entry, CHR_NUM_hg38=chrmsm, CHR_START_POS_hg38=startp, CHR_END_POS_hg38=endp, 
            CHR_START_POS_CODING_hg38=startpp, CHR_END_POS_CODING_hg38=endpp, CHR_STAND_hg38=strand_dr, 
            NT1_hg38=nt1,NT2_hg38=nt2,NT3_hg38=nt3)

        temp_Data <- Data[,-c(1,2)]
        
        all_Data <- cbind(tempf, temp_Data)
        fwrite(all_Data, paste0(store_dir,'/',basename(all_files[k])), row.names=FALSE, sep='\t', quote=FALSE)

    }
    
    cat('File',k,' of ', length(all_files), ' done\n')

}


