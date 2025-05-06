##############################################################################################
## Purpose: clinically relevant GRNs
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)
library(dorothea)
library(pheatmap)

save_dir <- '../results'
psi_input <- '../data/PSI_data'

fdr <- 0.05
diff <- 0.1
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
tf_ensemb_map <- as.data.frame(data.table::fread('../data/TF_ensembl_uniprot.txt', sep='\t'))

all_filesxx <- gtools::mixedsort(list.files(psi_input, pattern='*filtered_PSI.txt', full.names=TRUE))
all_cancer <- substr(basename(all_filesxx), 1,4)


##-- contradictory event pairs ------------------------------
cont_pairs <- data.table::fread('../data/contradictory_event_pair_DBD.txt')

##-- Contradictory event pairs that also show opposite survival ---- 
surv_pairs <- data.table::fread('../data/contradictory_event_survival.txt')

##-- Contradictory event-gene pairs ---- 
corr_pairs <- data.table::fread('../data/contradictory_event_pair-gene_comb.txt')

##--- Contradictory pairs showing contradictory survival associations ----
fdatay <- surv_pairs
num_sig_cors_d <- list()
for(k in 1:length(all_cancer)){
    fdatayt <- fdatay[fdatay$CANCER == all_cancer[k], ]
    fdatayu <- cont_pairs[cont_pairs$CANCER == all_cancer[k], ]
    fdatayt_gr <- igraph::graph_from_data_frame(fdatayt[,c(2,3)], directed=FALSE)
    fdatayu_gr <- igraph::graph_from_data_frame(fdatayu[,c(2,3)], directed=FALSE)
    num_sig_cors_d[[k]] <- igraph::intersection(fdatayt_gr, fdatayu_gr)
}

##--- splicing perturbed survival associated regulatory interactions ---
fdatay <- corr_pairs
fdatay$X1 <- paste0(fdatay$TF, '_', fdatay$ASID1)
fdatay$X2 <- paste0(fdatay$TF, '_', fdatay$ASID2)
num_sig_cors_dx <- list()
for(k in 1:length(all_cancer)){
    fdatayt <- fdatay[fdatay$CANCER == all_cancer[k], ]
    fdatayu <- igraph::as_data_frame(num_sig_cors_d[[k]])
    fdatayu$X1 <- unlist(lapply(fdatayu$from, function(x) substr(x, 1, nchar(x)-3)))
	fdatayu$X2 <- unlist(lapply(fdatayu$to, function(x) substr(x, 1, nchar(x)-3)))
    fdatayt_gr <- igraph::graph_from_data_frame(fdatayt[,c(8,9)], directed=FALSE)
    fdatayu_gr <- igraph::graph_from_data_frame(fdatayu[,c(3,4)], directed=FALSE)
    num_sig_cors_dx[[k]] <- igraph::intersection(fdatayt_gr, fdatayu_gr)
}

##--- Number of TF-gene pairs at least two splicing events of TF are contradictory in terms of 
## (1) splicing, (2)  survival, and (3) correlations with expression of the gene ---
num_sig_cors_dy <- list()
num_sig_cors_dy_u <- list()

for(k in 1:length(all_cancer)){
	
	temp <- igraph::as_data_frame(num_sig_cors_dx[[k]])

	corr_pairsc <- corr_pairs[corr_pairs$CANCER == all_cancer[k], ]
	corr_pairsc$X1 <- paste0(corr_pairsc$TF, '_', corr_pairsc$ASID1)
	corr_pairsc$X2 <- paste0(corr_pairsc$TF, '_', corr_pairsc$ASID2)

	pos <- c()
    for(j in 1:length(temp[[1]])){
        wh1 <- which(corr_pairsc$X1 == temp[[1]][j])
        wh2 <- which(corr_pairsc$X2 == temp[[2]][j])
        wha <- intersect(wh1, wh2)
        wh1 <- which(corr_pairsc$X2 == temp[[1]][j])
        wh2 <- which(corr_pairsc$X1 == temp[[2]][j])
        whb <- intersect(wh1, wh2)
        wh <- union(wha, whb)
        pos <- c(pos, wh)
    }

	corr_pairscx <- corr_pairsc[pos,]

	cont_pairsc <- cont_pairs[cont_pairs$CANCER == all_cancer[k], ]
	cont_pairsc$X1 <- unlist(lapply(cont_pairsc$AS1, function(x) substr(x, 1, nchar(x)-3)))
	cont_pairsc$X2 <- unlist(lapply(cont_pairsc$AS2, function(x) substr(x, 1, nchar(x)-3)))

	diff1 <- c()
	diff2 <- c()
	for(j in 1:length(corr_pairscx[[1]])){
		wh1 <- which(cont_pairsc$X1 == corr_pairscx$X1[j])
		wh2 <- which(cont_pairsc$X2 == corr_pairscx$X2[j])
		wha <- intersect(wh1, wh2)

		wh1 <- which(cont_pairsc$X1 == corr_pairscx$X2[j])
		wh2 <- which(cont_pairsc$X2 == corr_pairscx$X1[j])
		whb <- intersect(wh1, wh2)

		if(length(wha) !=0 ){
			diff1 <- c(diff1, cont_pairsc$MD1[wha])
			diff2 <- c(diff2, cont_pairsc$MD2[wha])
		}

		if(length(whb) !=0 ){
			diff1 <- c(diff1, cont_pairsc$MD2[whb])
			diff2 <- c(diff2, cont_pairsc$MD1[whb])
		}

	}

	surv_pairsc <- surv_pairs[surv_pairs$CANCER == all_cancer[k], ]
	surv_pairsc$X1 <- unlist(lapply(surv_pairsc$AS1, function(x) substr(x, 1, nchar(x)-3)))
	surv_pairsc$X2 <- unlist(lapply(surv_pairsc$AS2, function(x) substr(x, 1, nchar(x)-3)))

	hr1 <- c()
	hr2 <- c()
	for(j in 1:length(corr_pairscx[[1]])){

		wh1 <- which(surv_pairsc$X1 == corr_pairscx$X1[j])
		wh2 <- which(surv_pairsc$X2 == corr_pairscx$X2[j])
		wha <- intersect(wh1, wh2)

		wh1 <- which(surv_pairsc$X1 == corr_pairscx$X2[j])
		wh2 <- which(surv_pairsc$X2 == corr_pairscx$X1[j])
		whb <- intersect(wh1, wh2)

		if(length(wha) !=0 ){
			hr1 <- c(hr1, surv_pairsc$MD1[wha])
			hr2 <- c(hr2, surv_pairsc$MD2[wha])
		}

		if(length(whb) !=0 ){
			hr1 <- c(hr1, surv_pairsc$MD2[whb])
			hr2 <- c(hr2, surv_pairsc$MD1[whb])
		}

		if(length(wha) ==0 & length(whb) ==0){print(j)}

	}

	corr_pairscx$DIFF1 <- diff1
	corr_pairscx$DIFF2 <- diff2
	corr_pairscx$HR1 <- hr1
	corr_pairscx$HR2 <- hr2

	num_sig_cors_dy[[k]] <- corr_pairscx
	num_sig_cors_dy_u[[k]] <- igraph::simplify(igraph::graph_from_data_frame(corr_pairscx[,c(2,3)]))
}

pdata <- data.frame(x=all_cancer, freq=unlist(lapply(num_sig_cors_dy_u, function(x) igraph::ecount(x))))
p <- ggplot(pdata, aes(x, freq)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of TF-gene pairs with at least\n one contradictory event pair with \ncontradictory survival associations", limits=c(0,max(pdata$freq)+10)) +
geom_text(aes(label=freq), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0("../results/Survival_TF_gene_pairs_with_contradictory_events.png"),width=3.5, height=3, dpi=400)





##--- Contradictory pairs ----
cont_pairs_gr <- list()
for(k in 1:length(all_cancer)){
    fdatayu <- cont_pairs[cont_pairs$CANCER == all_cancer[k], ]
    cont_pairs_gr[[k]] <- igraph::graph_from_data_frame(fdatayu[,c(2,3)], directed=FALSE)
}

##--- splicing perturbed regulatory interactions ---
fdatay <- corr_pairs
fdatay$X1 <- paste0(fdatay$TF, '_', fdatay$ASID1)
fdatay$X2 <- paste0(fdatay$TF, '_', fdatay$ASID2)
cont_pairs_gr_dx <- list()
for(k in 1:length(all_cancer)){
    fdatayt <- fdatay[fdatay$CANCER == all_cancer[k], ]
    fdatayu <- igraph::as_data_frame(cont_pairs_gr[[k]])
    fdatayu$X1 <- unlist(lapply(fdatayu$from, function(x) substr(x, 1, nchar(x)-3)))
	fdatayu$X2 <- unlist(lapply(fdatayu$to, function(x) substr(x, 1, nchar(x)-3)))
    fdatayt_gr <- igraph::graph_from_data_frame(fdatayt[,c(8,9)], directed=FALSE)
    fdatayu_gr <- igraph::graph_from_data_frame(fdatayu[,c(3,4)], directed=FALSE)
    cont_pairs_gr_dx[[k]] <- igraph::intersection(fdatayt_gr, fdatayu_gr)
}

##--- Number of TF-gene pairs at least two splicing events of TF are contradictory in terms of 
## (1) splicing and (2) correlations with expression of the gene ---
cont_pairs_gr_dy <- list()
cont_pairs_gr_dy_u <- list()
for(k in 1:length(all_cancer)){
	temp <- igraph::as_data_frame(cont_pairs_gr_dx[[k]])
	corr_pairsc <- corr_pairs[corr_pairs$CANCER == all_cancer[k], ]
	corr_pairsc$X1 <- paste0(corr_pairsc$TF, '_', corr_pairsc$ASID1)
	corr_pairsc$X2 <- paste0(corr_pairsc$TF, '_', corr_pairsc$ASID2)

	pos <- c()
    for(j in 1:length(temp[[1]])){
        wh1 <- which(corr_pairsc$X1 == temp[[1]][j])
        wh2 <- which(corr_pairsc$X2 == temp[[2]][j])
        wha <- intersect(wh1, wh2)
        wh1 <- which(corr_pairsc$X2 == temp[[1]][j])
        wh2 <- which(corr_pairsc$X1 == temp[[2]][j])
        whb <- intersect(wh1, wh2)
        wh <- union(wha, whb)
        pos <- c(pos, wh)
    }

	corr_pairscx <- corr_pairsc[pos,]

	cont_pairsc <- cont_pairs[cont_pairs$CANCER == all_cancer[k], ]
	cont_pairsc$X1 <- unlist(lapply(cont_pairsc$AS1, function(x) substr(x, 1, nchar(x)-3)))
	cont_pairsc$X2 <- unlist(lapply(cont_pairsc$AS2, function(x) substr(x, 1, nchar(x)-3)))

	diff1 <- c()
	diff2 <- c()
	for(j in 1:length(corr_pairscx[[1]])){
		wh1 <- which(cont_pairsc$X1 == corr_pairscx$X1[j])
		wh2 <- which(cont_pairsc$X2 == corr_pairscx$X2[j])
		wha <- intersect(wh1, wh2)

		wh1 <- which(cont_pairsc$X1 == corr_pairscx$X2[j])
		wh2 <- which(cont_pairsc$X2 == corr_pairscx$X1[j])
		whb <- intersect(wh1, wh2)

		if(length(wha) !=0 ){
			diff1 <- c(diff1, cont_pairsc$MD1[wha])
			diff2 <- c(diff2, cont_pairsc$MD2[wha])
		}

		if(length(whb) !=0 ){
			diff1 <- c(diff1, cont_pairsc$MD2[whb])
			diff2 <- c(diff2, cont_pairsc$MD1[whb])
		}

	}

	corr_pairscx$DIFF1 <- diff1
	corr_pairscx$DIFF2 <- diff2

	cont_pairs_gr_dy[[k]] <- corr_pairscx
	cont_pairs_gr_dy_u[[k]] <- igraph::simplify(igraph::graph_from_data_frame(corr_pairscx[,c(2,3)]))
}

pdata <- data.frame(x=all_cancer, freq=unlist(lapply(cont_pairs_gr_dy_u, function(x) igraph::ecount(x))))
p <- ggplot(pdata, aes(x, freq)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of TF-gene pairs with at least\n one contradictory event pair", limits=c(0,max(pdata$freq)+50)) +
geom_text(aes(label=freq), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0("../results/TF_gene_pairs_with_contradictory_events.png"),width=3.5, height=3, dpi=400)

