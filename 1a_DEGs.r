##############################################################################################
# Purpose: Detect differentially expressed genes in each cancer type considering samples
# present in both the TCGA and splice data
# Should be adjusted to use hg19 build --> currently its hg38  --> Should it be?
##############################################################################################

rm(list=ls())
library(GenomicDataCommons)
library(DESeq2)
library(scales)
library(viridis)
library(ggplot2)
library(data.table)
library(biomaRt)

##---- Parameters --------
count <- 10 # number of patients with paired cancer and normal samples
##-----------------

TCGAtranslateID <- function(file_ids) {
    info <- GenomicDataCommons::files() %>%
        GenomicDataCommons::filter( ~ file_id %in% file_ids) %>%
        GenomicDataCommons::select('cases.samples.submitter_id') %>%
        results_all()

    # The mess of code below is to extract TCGA barcodes
    # id_list will contain a list (one item for each file_id)
    # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
    id_list <- lapply(info$cases,function(a) {
        a[[1]][[1]][[1]]})
    # so we can later expand to a data.frame of the right size
    barcodes_per_file <- sapply(id_list,length)
    # And build the data.frame
    return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                      submitter_id = unlist(id_list)))
}

TCGAconsideredCancerTypes <- function(){

    ##----- Number of cases per project ------------------------------------
    res <- GenomicDataCommons::cases() %>% GenomicDataCommons::facet("project.project_id") %>% GenomicDataCommons::aggregations()
    res1 <- res$project.project_id
    res1$temp <- substr(res1[[2]], 1,4)
    res2 <- res1[res1$temp == 'TCGA', ]
    temp_type <- unlist(lapply(strsplit(res2$key, '[-]'),'[[',2))

    ##---- to store the cancer-types with at least 10 samples with both normal and cancer tissue
    c_type <- c()

    ##--- check for number of samples in both cancer and normal --------
    for(k in 1:length(temp_type)){

        ##----- Get manifest -----------------------------------------------
        ge_manifest = GenomicDataCommons::files() %>%
            GenomicDataCommons::filter( cases.project.project_id == paste0('TCGA-',temp_type[k])) %>% 
            GenomicDataCommons::filter( type == 'aligned_reads' ) %>%
            GenomicDataCommons::manifest()
        ge_manifest1 <- as.data.frame(ge_manifest)

        ##---- only RNA seq aligned to genome data
        temp_bam <- ge_manifest1
        wh3 <- which(temp_bam$file_name %like% 'rna_seq.genomic')
        temp_bam1 <- temp_bam[wh3, ]

        temp_bam1 <- temp_bam1[,c(1,19,12)]
        colnames(temp_bam1) <- c('file_id','filename','md5')
        sids <- TCGAtranslateID(temp_bam1[[1]])
        temp_bam2 <- merge(temp_bam1, sids, by='file_id')

        temp_bam2$nid <- substr(temp_bam2$submitter_id, 1, 12)
        # only keep normal tissue (11) and primary cancer (01)
        temp_bam2$sample_id <- substr(temp_bam2$submitter_id, 14,15)
        temp_bam3 <- temp_bam2[temp_bam2$sample_id %in% c('01','11'), ]

        ## choose only those submitter ids that have at least one in each category 01 and 11
        tempx1 <- temp_bam3[temp_bam3$sample_id == '01', ]$nid
        tempx2 <- temp_bam3[temp_bam3$sample_id == '11', ]$nid
        tempx3 <- intersect(tempx1, tempx2)
        temp_bam4 <- temp_bam3[temp_bam3$nid %in% tempx3, ]

        # only keep the patients with one sample in each of normal solid tissue and primary cancer
        allnids <- unique(temp_bam4$nid)
        whc <- c()
        for(j in 1:length(allnids)){
            temp_wh <- which(temp_bam4$nid == allnids[j])
            temp <- temp_bam4[temp_wh, ]
            if(length(temp_wh) == 2){
                whc <- c(whc, temp_wh)
            }else{
                temp_wh1 <- which(temp$sample_id == '01')
                whc <- c(whc, temp_wh[temp_wh1[1]])

                temp_wh1 <- which(temp$sample_id == '11')
                whc <- c(whc, temp_wh[temp_wh1[1]])
            }
        }

        temp_final <- temp_bam4[whc, ]
        tcount <- count*2
        if(length(temp_final[[1]]) >= tcount){
            c_type <- c(c_type, temp_type[k])
        }
    }

    return(c_type)
}

c_select <- TCGAconsideredCancerTypes() ##-- 15 cancer types selected -----------
c_select <- gtools::mixedsort(c_select)

##---- Download expression counts of the samples --------------------------------------------
mani_dir <- '../data/Manifests'
if(!dir.exists(mani_dir)){
    dir.create(mani_dir)
}

expr_dir <- '../../../public_data/TCGA_expressions'
if(!dir.exists(expr_dir)){
    dir.create(expr_dir)
}

##---- downloaded files -------------------------
dwnld_dirs <- basename(list.files(expr_dir))
dwnld_files <- substr(basename(list.files(expr_dir, pattern='*.tsv.partial', recursive=TRUE)),1,36)

# counter <- 0
for(k in 1:length(c_select)){

    ge_manifest = files() %>%
        filter( cases.project.project_id == paste0('TCGA-',c_select[k])) %>% 
        filter( type == 'aligned_reads' ) %>%
        manifest()
    ge_manifest1 <- as.data.frame(ge_manifest)

    ge_manifest <- files() |>
        filter( cases.project.project_id == paste0('TCGA-',c_select[k])) |> 
        filter( type == 'gene_expression' ) |>
        filter( analysis.workflow_type == 'STAR - Counts')  |>
    manifest()
    ge_manifest1 <- as.data.frame(ge_manifest)
    ge_manifest2 <- ge_manifest1[ge_manifest1$data_type == 'Gene Expression Quantification',]
    ge_manifest3 <- ge_manifest2[ge_manifest2$access == 'open', ]
    ge_manifest3$file_name_st <- substr(ge_manifest3$file_name,1,36)

    temp_file <- paste0(mani_dir,'/',c_select[k],'_manifest_final.txt')
    data.table::fwrite(ge_manifest3, temp_file, sep='\t', row.names=FALSE, quote=FALSE)

    # counter <- counter + length(ge_manifest3$file_id)
    whs <- setdiff(ge_manifest3$file_id, dwnld_dirs)

    if(length(whs) != 0){
        tempx <- ge_manifest3[ge_manifest3$file_id %in% whs, ]
        temp_filex <- paste0('../data/manifest_final_temp.txt')
        data.table::fwrite(tempx, temp_filex, sep='\t', row.names=FALSE, quote=FALSE)
        temp_cmd <- paste0('./gdc-client download -n 5 -d ',expr_dir,' -m ', temp_filex)
        system(temp_cmd)
    }

    whs <- intersect(ge_manifest3$file_name_st, dwnld_files)

    if(length(whs) != 0){
        tempx <- ge_manifest3[ge_manifest3$file_name_st %in% whs, ]
        temp_filex <- paste0('../data/manifest_final_temp.txt')
        data.table::fwrite(tempx, temp_filex, sep='\t', row.names=FALSE, quote=FALSE)
        temp_cmd <- paste0('./gdc-client download -n 5 -d ',expr_dir,' -m ', temp_filex)
        system(temp_cmd)
    }
    cat('Cancer',k,'of',length(c_select),'done\n')
}

###------------------------------------------------------------------------------------------------




#####------- Write the gene counts in proper format ----------------------------------------------

mani_dir <- '../data/Manifests'
manifests <- list.files(mani_dir, pattern='*_final.txt', full.names=TRUE)
all_psi <- list.files('../data/PSI_data', pattern='^PSI_download', full.names=TRUE)
expr_dir <- '../../../public_data/TCGA_expressions'
out_dir <- '../../../public_data/TCGA_expr_counts'
dir.create(out_dir)

for(k in 1:length(manifests)){

    temp <- data.table::fread(manifests[k])
    temp1 <- temp[,c(1,4,12)]
    colnames(temp1) <- c('file_id','filename','md5')
    sids <- TCGAtranslateID(temp1[[1]])
    temp2 <- merge(temp1, sids, by='file_id')
    temp2$sample_num <- substr(temp2$submitter_id, 16,16)

    ## some patients have more than one normal or cancer samples. Selecting only patients with unique cancer and normal samples
    ## This is because of the following reasons.
    ## In the TCGA splice site only one normal/cancer sample per patient is included. In this sense, if a patient had multiple normal
    ## samples, for example, then I do not know which of the sample was used for the splicing calculations
    temp2 <- temp2[temp2$sample_num == 'A', ] 
    temp2$submitter_id <- substr(temp2$submitter_id, 1,15)
    templt <- plyr::count(temp2$submitter_id)
    lt <- templt[templt$freq > 1, ]$x 
    all_submitter <- temp2$submitter_id
    temp2 <- temp2[temp2$submitter_id %in% setdiff(all_submitter, lt), ]

    ##-- read the corresponding PSI file ------------------------------------------
    temp_cancer <- strsplit(basename(manifests[k]),'[_]')[[1]][1]
    wh <- which(all_psi %like% temp_cancer)
    tempf <- data.table::fread(all_psi[wh], sep='\t', fill=TRUE)
    len <- length(tempf) - 2
    temp_cols <- colnames(tempf)[11:len]
    temp_cols <- gsub('_','-',temp_cols)
    temp_cols <- gsub('Norm','11',temp_cols)
    wh <- which(nchar(temp_cols) == 12)
    temp_cols[wh] <- paste0(temp_cols[wh],'-01')

    ##---- select from the TCGA data ---
    temp3 <- temp2[temp2$submitter_id %in% temp_cols, ]

    temp_file <- paste0(mani_dir,'/',temp_cancer,'_manifest_PSI_and_Expr.txt')
    data.table::fwrite(temp3, temp_file, sep='\t', row.names=FALSE, quote=FALSE)

    tempy <- as.data.frame(data.table::fread(paste0(expr_dir,'/',temp3$file_id[1],'/',temp3$filename[1])))
    tempy <- tempy[5:length(tempy[[1]]),c(1,4)]
    tempz <- tempy[order(tempy$gene_id),]

    for(j in 2:length(temp3[[1]])){
        tempy <- as.data.frame(data.table::fread(paste0(expr_dir,'/',temp3$file_id[j],'/',temp3$filename[j])))
        tempy <- tempy[5:length(tempy[[1]]),c(1,4)] ## taking the Ensembl ID and unstranded counts
        tempyy <- tempy[order(tempy$gene_id),]
        tempz <- cbind(tempz, tempyy$unstranded)
    }

    colnames(tempz) <- c('Gene_Symbol',temp3$submitter_id)
    data.table::fwrite(tempz, paste0(out_dir,'/',temp_cancer,'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    cat('Cancer', k, 'of',length(manifests),'done\n')

}

#####-------------------------------------------------------------------------------------------------------------------

exons_f <- data.table::fread('../data/Ensembl_exons.txt', header=TRUE)
exons_f$len <- abs(exons_f$V5 - exons_f$V4)+1
all_genes <- unique(exons_f$gene_id)
lens <- c()
for(k in 1:length(all_genes)){
    lens <- c(lens, sum(exons_f[exons_f$gene_id == all_genes[k], ]$len))
    cat('Gene',k, 'of',length(all_genes),'done\n')
}
exlen <- data.frame(Gene_Symbol=all_genes, Length=lens)

##--------- DIfferential Gene Expression -------------------------------------------------------------------------------
gene_counts <- list.files(out_dir,pattern='*.txt', full.names=TRUE)
save_dir <- '../results/Diff_expr'
if(!dir.exists(save_dir)){
    dir.create(save_dir)
}

save_dirx <- '../data/Diff_expr'
if(!dir.exists(save_dirx)){
    dir.create(save_dirx)
}

save_diry <- '../../../public_data/TCGA_normalized_counts'
if(!dir.exists(save_diry)){
    dir.create(save_diry)
}

save_dirf <- '../../../public_data/TCGA_FPKM_counts'
if(!dir.exists(save_dirf)){
    dir.create(save_dirf)
}

for(k in 1:length(manifests)){

    temp_file <- as.data.frame(data.table::fread(gene_counts[k]))
    tgenes <- unlist(lapply(strsplit(temp_file[[1]],'[.]'),'[[',1))
    wh <- which(tgenes %in% all_genes)
    tgenes <- tgenes[wh]
    temp_file1 <- temp_file[wh, ]
    gnc <- plyr::count(tgenes)
    wh <- which(gnc$freq > 1)
    g2c <- setdiff(tgenes, gnc$x[wh])
    wh <- which(tgenes %in% g2c)
    tgenes <- tgenes[wh]
    temp_file1 <- temp_file[wh, ]
    temp_file1$Gene_Symbol <- tgenes
    temp_genes <- temp_file1[[1]]
    temp_file1 <- temp_file1[,-1]
    rownames(temp_file1) <- temp_genes
    temp_file <- temp_file1

    last_flag <- substr(colnames(temp_file), 14,15)
    n_col1 <- which(last_flag == '11') ## control samples
    n_col2 <- which(last_flag == '01') ## cancer samples

    ## remove gene names columns
    group <- rep("", length(last_flag))
    group[n_col1] <- 'Control'
    group[n_col2] <- 'Condition'

    ##-- prepare the sample data
    sample_info <- data.frame(group=group)
    rownames(sample_info) <- colnames(temp_file)

    xx <- all(colnames(temp_file) %in% rownames(sample_info)) ## all row and colnames present
    yy <- all(colnames(temp_file) == rownames(sample_info)) ## all row and colnames present in the same order

    if(!xx | !yy){
        break
    }

    ##--- DESeq2 format ---
    x.deseq2 <- DESeq2::DESeqDataSetFromMatrix(countData=temp_file, colData=sample_info, design= ~ group)

    ##--- Prefiltering ------------
    # dim(x.deseq2)
    x.deseq2.filt <- x.deseq2#[rowSums(DESeq2::counts(x.deseq2)) > 10, ] ## At least 10 reads in all samples
    # dim(x.deseq2.filt)
    ##-----------------------------
    ##-- set the reference to control ----
    x.deseq2.filt$group <- relevel(x.deseq2.filt$group, ref='Control')
    ##-----------------------------------
    
    ##---- normalized counts -----
    x.deseq2.filt.normalized <- DESeq2::estimateSizeFactors(x.deseq2.filt)
    x.deseq2.filt.normalized.counts <- as.data.frame(counts(x.deseq2.filt.normalized, normalized=TRUE))
    x.deseq2.filt.normalized.counts$Ensembl_gene_id <- rownames(x.deseq2.filt.normalized.counts)

    wh <- which(exlen$Gene_Symbol %in% x.deseq2.filt.normalized.counts$Ensembl_gene_id)
    exlenx <- exlen[wh, ]
    exlenx <- exlenx[match(x.deseq2.filt.normalized.counts$Ensembl_gene_id, exlenx$Gene_Symbol),]
    mcols(x.deseq2.filt.normalized)$basepairs <- exlenx$Length
    fpkm_data <- fpkm(x.deseq2.filt.normalized)
    # x.deseq2.filt.fpkm.counts <- as.data.frame(counts(fpkm_data, normalized=TRUE))
    # x.deseq2.filt.fpkm.counts$Ensembl_gene_id <- rownames(x.deseq2.filt.fpkm.counts)

    data.table::fwrite(x.deseq2.filt.normalized.counts, paste0(save_diry,"/",substr(basename(manifests[k]),1,4),"_normalized_counts.txt"), sep='\t', row.names=FALSE, quote=FALSE)
    data.table::fwrite(fpkm_data, paste0(save_dirf,"/",substr(basename(manifests[k]),1,4),"_fpkm_counts.txt"), sep='\t', row.names=FALSE, quote=FALSE)

    ##----------------------------
    x.deseq2.filt.normalized.de <- DESeq2::DESeq(x.deseq2.filt.normalized)

    res <- DESeq2::results(x.deseq2.filt.normalized.de)

    ##----- save and plot the results ------------------------
    deseq2ResDF <- as.data.frame(res)
    deseq2ResDF$Ensembl_gene_id <- rownames(deseq2ResDF)
    data.table::fwrite(deseq2ResDF, paste0(save_dirx,"/",substr(basename(manifests[k]),1,4),"_diff_expr.txt"), sep='\t', row.names=FALSE, quote=FALSE)

    # Set a boolean column for significance
    deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < 0.05, "Significant", NA)

    # Plot the results 
    pp <- ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + 
    geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + 
    scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + 
    labs(x="mean of normalized counts", y="log fold change") + 
    scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + 
    geom_density_2d(colour="black", size=1)

    ggsave(pp,filename=paste0(save_dir,"/",substr(basename(manifests[k]),1,4),"_diff_expr.png"),width=7, height=4, dpi=400)

    cat('Cancer',k,'of', length(manifests),'done\n')

}


