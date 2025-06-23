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

paired_sam <- c()
tcancer <- c()

##-- get patient IDs with paired samples-------
for(k in 1:length(manifests)){

    temp_mani <- data.table::fread(manifests[k])

    ##-- which patient id has both normal and cancer samples --
    pids <- plyr::count(substr(temp_mani$submitter_id, 1,12))
    pids_filt <- pids[which(pids$freq > 1),]$x
    pids_filt <- gsub('-','_',pids_filt)
    ##--

	##-- read the corresponding PSI file --------------------------
    c_select <- substr(basename(manifests[k]),1,4)
	wh <- which(all_psi %like% c_select)
	temp_psi <- data.table::fread(all_psi[wh], sep='\t', fill=TRUE)
	temp_psi <- as.data.frame(temp_psi)
	##-- patient IDs start from column 11, and splciing events start from row 14!

	col_pos <- c()
	for(j in 1:length(pids_filt)){
		wh <- which(colnames(temp_psi) %like% pids_filt[j])
		col_pos <- union(col_pos,wh)
	}

	##---- Some TCGA patients with paired samples are missing paired samples in the PSI data ---
	##---- detect which patient has no paired samples in the PSI data and remove it from the col_pos variable ---
	patient_ids <- substr(colnames(temp_psi)[col_pos], 1,12)
	patient_sample_counts <- plyr::count(patient_ids)
	to_remove_patient <- patient_sample_counts[[1]][which(patient_sample_counts$freq < 2)]
	whc <- which(patient_ids %in% to_remove_patient)
	if(length(whc) > 0){col_pos <- col_pos[-whc]}
	##------------------------------------------------------------
	paired_sam <- c(paired_sam, length(col_pos)/2)
	tcancer <- c(tcancer, c_select)

	##-- samples ----
	keep_cols <- c(seq(1,10),col_pos)
	temp_psi1 <- temp_psi[ ,keep_cols] ## contains first 10 columns and all sample columns but removes some columns from end which are not samples
	temp_psi2 <- temp_psi[14:length(temp_psi[[1]]), col_pos] # only keeps the samples columns and splicing event rows

	###--- divide into normal and cancer samples ---
	normal_id <- which(colnames(temp_psi2) %like% '_Norm')
	cancer_id <- setdiff(seq(1,length(temp_psi2)), normal_id)
	normal <- temp_psi2[, normal_id]
	cancer <- temp_psi2[, cancer_id]
	normal <- normal[,order(colnames(normal))]
	cancer <- cancer[,order(colnames(cancer))]

	wh1 <- apply(normal, 1, function(x) which(!is.na(as.numeric(x))))
	wh2 <- apply(cancer, 1, function(x) which(!is.na(as.numeric(x))))
	wh3 <- mapply(function(x,y) intersect(x,y), wh1, wh2)
	present_both <- which(lengths(wh3)  >= count)

	normal <- normal[present_both, ]
	cancer <- cancer[present_both, ]
	whx <- wh3[present_both]

	normal_list <- as.list(as.data.frame(t(normal)))
	cancer_list <- as.list(as.data.frame(t(cancer)))

	normal_list_filt <- mapply(function(x,y) as.numeric(x[y]), normal_list, whx)
	cancer_list_filt <- mapply(function(x,y) as.numeric(x[y]), cancer_list, whx)

	can_norm_diff <- mapply(function(x,y) as.numeric(x) - as.numeric(y), cancer_list_filt, normal_list_filt)
	allen <- lengths(can_norm_diff)
	pos_len <- lapply(can_norm_diff, function(x) length(which(x >= 0)))
	neg_len <- lapply(can_norm_diff, function(x) length(which(x < 0)))
	# pos_flag <- mapply(function(x,y) x/y, pos_len, allen)
	# pos_wh <- which(pos_flag == 1)
	# neg_flag <- mapply(function(x,y) x/y, neg_len, allen)
	# neg_wh <- which(neg_flag == 1)


	pvals <- mapply(function(x,y) wilcox.test(as.numeric(x), as.numeric(y), paired=TRUE)$p.value, normal_list_filt, cancer_list_filt)
	
	diff_normal_cancer <- mapply(function(x,y) as.numeric(x)-as.numeric(y), cancer_list_filt, normal_list_filt)
	median_normal <- lapply(normal_list_filt, function(x) median(as.numeric(x)))
	mean_normal <- lapply(normal_list_filt, function(x) mean(as.numeric(x)))
	median_cancer <- lapply(cancer_list_filt, function(x) median(as.numeric(x)))
	mean_cancer <- lapply(cancer_list_filt, function(x) mean(as.numeric(x)))
	diff_median <- lapply(diff_normal_cancer, function(x) median(as.numeric(x)))
	diff_mean <- lapply(diff_normal_cancer, function(x) mean(as.numeric(x)))
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
	temp_psi3$POS <- pos_len
	temp_psi3$NEG <- neg_len

	wh2 <- present_both+13 ## because we need to select these rows from temp_psi1
	keep_rows <- c(seq(1:13),wh2)
	temp_psi1 <- temp_psi1[keep_rows,]

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
	tempx$POS <- rep(NA, length(tempx[[1]]))
	tempx$NEG <- rep(NA, length(tempx[[1]]))
	colnames(tempx) <- colnames(tempy)
	temp_psi5 <- rbind(tempx, tempy)

	data.table::fwrite(temp_psi5, paste0('../data/PSI_data/',c_select,'_filtered_PSI_paired.txt'), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
	cat('Cancer',k,'of',length(manifests),'done\n')

}

data.table::fwrite(data.frame(CANCER=tcancer, PAIRED_SAMPLES=paired_sam), paste0('../data/cancer_paired_samples.txt'), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)


