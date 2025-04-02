##############################################################################################
# Purpose: Save the filtered PSI value files --
##############################################################################################

rm(list=ls())
library(data.table)
# library(GenomicDataCommons)

##---- Parameters --------
count <- 10 # number of patients in cancer and normal samples for which a splicing event should have a valid value
##-----------------

mani_dir <- '../data/Manifests' ## manifests filtered as per TCGA gene expression availability

all_psi <- list.files('../data/PSI_data', pattern='^PSI_download', full.names=TRUE)
manifests <- gtools::mixedsort(list.files(mani_dir, pattern = '*_manifest_PSI_and_Expr.txt', full.names=TRUE))


##-- get patient IDs with paired samples-------
for(k in 1:length(manifests)){

    temp_mani <- data.table::fread(manifests[k])

	##-- read the corresponding PSI file --------------------------
    c_select <- substr(basename(manifests[k]),1,4)
	wh <- which(all_psi %like% c_select)
	temp_psi <- data.table::fread(all_psi[wh], sep='\t', fill=TRUE)
	temp_psi <- as.data.frame(temp_psi)
	##-- patient IDs start from column 11, and splciing events start from row 14!

	##-- samples ----
	col_pos <- which(colnames(temp_psi) %like% 'TCGA')
	keep_cols <- c(seq(1,10),col_pos)
	temp_psi1 <- temp_psi[ ,keep_cols] ## contains first 10 columns and all sample columns but removes some columns from end which are not samples

	temp_psi2 <- temp_psi[14:length(temp_psi[[1]]), col_pos] # only keeps the samples columns and splicing event rows

    ##-- Filter to only retain samples also present in the TCGA gene expression data
    temp_mani_sid <- gsub('-','_',substr(temp_mani$submitter_id, 1, 12))
    temp_psi2_sid <- substr(colnames(temp_psi2), 1, 12)
    wh_samples <- which(temp_psi2_sid %in% temp_mani_sid)
    temp_psi2 <- temp_psi2[,wh_samples]


	###--- divide into normal and cancer samples ---
	normal_id <- which(colnames(temp_psi2) %like% '_Norm')
	cancer_id <- setdiff(seq(1,length(temp_psi2)), normal_id)
	normal <- temp_psi2[, normal_id]
	cancer <- temp_psi2[, cancer_id]

	##-- only keep the splicing events that have values for at least count (here, 10) number of cancer and normal samples --
	wh <- apply(normal, 1, function(x) length(which(x != "null")))
	present_normal <- which(wh  >= count)
	wh <- apply(cancer, 1, function(x) length(which(x != "null")))
	present_cancer <- which(wh  >= count)
	present_both <- intersect(present_cancer, present_normal)
	normal <- normal[present_both, ]
	cancer <- cancer[present_both, ]

	normal_list <- as.list(as.data.frame(t(normal[,order(colnames(normal))])))
	normal_list_filt <- lapply(normal_list, function(x) x[x != 'null']) ## for each event contains the value for the normal samples
	cancer_list <- as.list(as.data.frame(t(cancer[,order(colnames(cancer))])))
	cancer_list_filt <- lapply(cancer_list, function(x) x[x != 'null']) ## for each event contains the value for the cancer samples

	pvals <- mapply(function(x,y) wilcox.test(as.numeric(x), as.numeric(y), paired=FALSE)$p.value, normal_list_filt, cancer_list_filt)
	median_normal <- lapply(normal_list_filt, function(x) median(as.numeric(x)))
	mean_normal <- lapply(normal_list_filt, function(x) mean(as.numeric(x)))
	median_cancer <- lapply(cancer_list_filt, function(x) median(as.numeric(x)))
	mean_cancer <- lapply(cancer_list_filt, function(x) mean(as.numeric(x)))
	diff_median <- mapply(function(x,y) median(as.numeric(y)) - median(as.numeric(x)), normal_list_filt, cancer_list_filt)
	diff_mean <- mapply(function(x,y) mean(as.numeric(y)) - mean(as.numeric(x)), normal_list_filt, cancer_list_filt)
	qvals <- p.adjust(pvals, 'fdr')

	##-- combine the normal and cancer data sets in one data frame
	temp_psi3 <- cbind(normal, cancer)
	temp_psi3$p_value <- pvals
	temp_psi3$FDR <- qvals
	temp_psi3$MEDIAN_DIFF <- diff_median
	temp_psi3$MEDIAN_NORMAL <- median_normal
	temp_psi3$MEDIAN_CANCER <- median_cancer
	temp_psi3$MEAN_DIFF <- diff_mean
	temp_psi3$MEAN_NORMAL <- mean_normal
	temp_psi3$MEAN_CANCER <- mean_cancer

	wh2 <- present_both+13 ## because we need to select these rows from temp_psi1
    wh3 <- wh_samples+10 ## because we need to select these samples from temp_psi1

	# # keep_rows <- c(seq(1:13),wh2)
    # # keep_cols <- c(seq(1:10),wh3)
	# temp_psi1a <- temp_psi1[1:13,1:10]
	# temp_psi1b <- temp_psi1[wh2,wh3]
	keep_rows <- c(seq(1:13),wh2)
    keep_cols <- c(seq(1:10),wh3)
	temp_psi1 <- temp_psi1[keep_rows,keep_cols]

	##----------------------------------------------------------
	tempy <- cbind(temp_psi1[14:length(temp_psi1[[1]]), 1:10], temp_psi3)
	tempx <- temp_psi1[1:13,]
	tempx$p_value <- rep(NA, length(tempx[[1]]))
	tempx$FDR <- rep(NA, length(tempx[[1]]))
	tempx$MEDIAN_DIFF <- rep(NA, length(tempx[[1]]))
	tempx$MEDIAN_NORMAL <- rep(NA, length(tempx[[1]]))
	tempx$MEDIAN_CANCER <- rep(NA, length(tempx[[1]]))
	tempx$MEAN_DIFF <- rep(NA, length(tempx[[1]]))
	tempx$MEAN_NORMAL <- rep(NA, length(tempx[[1]]))
	tempx$MEAN_CANCER <- rep(NA, length(tempx[[1]]))
	colnames(tempx) <- colnames(tempy)
	temp_psi5 <- rbind(tempx, tempy)

	data.table::fwrite(temp_psi5, paste0('../data/PSI_data/',c_select,'_filtered_PSI.txt'), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
	cat('Cancer',k,'of',length(manifests),'done\n')

}



