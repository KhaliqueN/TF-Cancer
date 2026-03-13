##############################################################################################
# Purpose: Save the filtered PSI value files --
##############################################################################################

rm(list=ls())
library(data.table)
#####---- From TCGA, we gathered that the following 15 cancer types having 10 matched cancer and normal samples with RNA-seq reads:
##### 'BLCA', 'BRCA', 'COAD','ESCA','HNSC','KICH',  'KIRC',  'KIRP','LIHC', 'LUAD', 'LUSC','PRAD','STAD','THCA', 'UCEC' ----
#### PSI files for these cancer types were manually downloaded from TCGASPliceSeq database ------------

save_dir <- 'results_rep/'
if(!dir.exists(save_dir)){
	dir.create(save_dir, recursive=TRUE)
}

count <- 10 ## least number of paired normal/cancer samples in which an event has valid values
all_psi <- list.files('data/PSI_data', pattern='^PSI_download', full.names=TRUE)
paired_sam <- c()
tcancer <- c()
asevents <- c()
can_samples <- c()

for(k in 1:length(all_psi)){

	##-- read the corresponding PSI file --------------------------
	temp_psi <- data.table::fread(all_psi[k], sep='\t', fill=TRUE)
	temp_psi <- as.data.frame(temp_psi)
	##-- patient IDs start from column 11, and splciing events start from row 14!

	temp_col_pos1 <- which(colnames(temp_psi) %like% '_Norm')
	temp_col_pos2 <- which(colnames(temp_psi) %like% 'TCGA_')
	temp_col_pos <- intersect(temp_col_pos1, temp_col_pos2)

	pids_filt <- unique(substr(colnames(temp_psi)[temp_col_pos], 1, 12))
	can_samples <- c(can_samples, length(which(colnames(temp_psi) %like% "TCGA_")) - length(temp_col_pos)) # store # of cancer samples

	col_pos <- c() ## collect all paired normal and cancer samples column positions in the data file 
	for(j in 1:length(pids_filt)){
		wh1 <- which(colnames(temp_psi) == paste0(pids_filt[j], '_Norm'))[1]
		wh2 <- which(colnames(temp_psi) == pids_filt[j])
		if(length(wh2) != 0){
			wh2 <- wh2[1]
		}else{
			next
		}
		col_pos <- union(col_pos, union(wh1, wh2))
	}

	paired_sam <- c(paired_sam, length(col_pos)/2) # store # of paired samples
	c_select <- substr(basename(all_psi[k]), 14, 17)
	tcancer <- c(tcancer, c_select) # stire cancer type

	##-- samples ----
	keep_cols <- c(seq(1,10), col_pos)
	temp_psi1 <- temp_psi[ ,keep_cols] ## contains first 10 columns and all sample columns but removes some columns from end which are not samples
	temp_psi2 <- temp_psi[14:length(temp_psi[[1]]), col_pos] # only keeps the samples columns and splicing event rows

	###--- divide into normal and cancer samples ---
	normal_id <- which(colnames(temp_psi2) %like% '_Norm')
	cancer_id <- setdiff(seq(1,length(temp_psi2)), normal_id)
	normal <- temp_psi2[, normal_id]
	cancer <- temp_psi2[, cancer_id]
	normal <- normal[,order(colnames(normal))]
	cancer <- cancer[,order(colnames(cancer))]

	## for each splicing event select samples with valid (not "null") values
	wh1 <- apply(normal, 1, function(x) which(!is.na(as.numeric(x))))
	wh2 <- apply(cancer, 1, function(x) which(!is.na(as.numeric(x))))
	wh3 <- mapply(function(x,y) intersect(x,y), wh1, wh2)
	present_both <- which(lengths(wh3)  >= count) ## only keeping the events with valid values for at least 10 paired cancer /normal samples

	normal <- normal[present_both, ]
	cancer <- cancer[present_both, ]
	whx <- wh3[present_both]

	normal_list <- as.list(as.data.frame(t(normal)))
	cancer_list <- as.list(as.data.frame(t(cancer)))

	normal_list_filt <- mapply(function(x,y) as.numeric(x[y]), normal_list, whx)
	cancer_list_filt <- mapply(function(x,y) as.numeric(x[y]), cancer_list, whx)

	can_norm_diff <- mapply(function(x,y) x - y, cancer_list_filt, normal_list_filt) 

	allen <- lengths(can_norm_diff)
	pos_len <- lapply(can_norm_diff, function(x) length(which(x >= 0)))
	neg_len <- lapply(can_norm_diff, function(x) length(which(x < 0)))

	## Here we perform wilxocon signed-rank test
	pvals <- mapply(function(x,y) wilcox.test(as.numeric(x), as.numeric(y), paired=TRUE)$p.value, normal_list_filt, cancer_list_filt)
	median_normal <- lapply(normal_list_filt, function(x) median(as.numeric(x)))
	mean_normal <- lapply(normal_list_filt, function(x) mean(as.numeric(x)))
	median_cancer <- lapply(cancer_list_filt, function(x) median(as.numeric(x)))
	mean_cancer <- lapply(cancer_list_filt, function(x) mean(as.numeric(x)))
	diff_median <- lapply(can_norm_diff, function(x) median(as.numeric(x)))
	diff_mean <- lapply(can_norm_diff, function(x) mean(as.numeric(x))) ## delta PSI
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
	asevents <- c(asevents, length(temp_psi5[[1]])-13) ## number of events with valid values for at keast 10 paired samples

	data.table::fwrite(temp_psi5, paste0('data/PSI_data/',c_select,'_PSI_paired.txt'), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
	cat('Cancer',k,'of',length(all_psi),'done\n')

}

tempd <- data.frame(CANCER=tcancer, PAIRED_SAMPLES=paired_sam, CANCER_SAMPLES=can_samples, AS_EVENTS=asevents)
data.table::fwrite(tempd, paste0(save_dir,'cancer_paired_samples.txt'), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

data.table::fwrite(tempd, paste0('data/cancer_paired_samples.txt'), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)


