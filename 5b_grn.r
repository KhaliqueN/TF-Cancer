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
cont_pairs <- data.table::fread(paste0('../data/contradictory_event_pair_DBD.txt'))

##-- Contradictory event pairs that also show opposite survival ---- 


##-- Contradictory event-gene pairs ---- 












# all_files_genes <- gtools::mixedsort(list.files('../../../public_data/TCGA_normalized_counts', full.names=TRUE))
# ## Some of the event coordinates from TCGA splice seq does not overlap to any of the exons mapped from the canonical protein
# ## sequence from Uniprot. This is because the canonical seqeunce is not always the longest.

# #####----- Using Dorothea for GRN --------------
# all_map <- data.table::fread('../data/ensembl_name_map.txt')
# all_mapx <- unique(all_map[,c(1,12,13)])
# conf_grn <- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c('A','B','C'), ])
# ##------------- 13,223 edges -------------------
# colnames(conf_grn) <- c('HGNC_symbol','Confidence','Target','Direction')
# conf_grnx <- merge(conf_grn, all_mapx, by='HGNC_symbol')
# conf_grnx <- conf_grnx[conf_grnx$Uniprotswissprot != '', ]
# colnames(conf_grnx) <- c('TF','Confidence','HGNC_symbol','Direction','Ensembl_gene_id','Uniprotswissprot')
# conf_grny <- merge(conf_grnx, all_mapx, by='HGNC_symbol')
# colnames(conf_grny) <- c('Target_symbol','TF_symbol','Confidence','Direction','TF_Ensembl_gene_id','TF_Uniprotswissprot',
#     'Target_Ensembl_gene_id','Target_Uniprotswissprot')
# conf_grny <- conf_grny[conf_grny$Target_Uniprotswissprot != '', ]
# conf_grnz <- conf_grny[complete.cases(conf_grny), ]
# ### 13,108 edges remain as per gene symbols --- 13,221 edges remain as per uniprot --- 15,243 edges remain as per Ensembl


# ##-- Taking gene symbol based ---
# conf_grn_symbol <- unique(conf_grnz[,c(1,2,3,4)])


# ##-- Take the edges concerning pancancer perturbed events ----
# tdatay <- readRDS('../data/Contradictory_event_pairs.rds')## not filtered by DBD perturbation
# # tdatay <- readRDS('../data/Contradictory_event_pairs_DBD.rds') ## filtered by DBD perturbation
# tempdata <- tdatay[[1]]  ## list of event pairs 
# tempx <- tdatay[[2]]  ## overlap of event pairs

# tedges <- data.frame(igraph::as_edgelist(tempx[[10]]))
# tnodes <- unique(unlist(lapply(strsplit(tedges[[1]],'[_]'),'[[',1))) ## names of the TFs perturbed in at least 5 cancer types
# tnet1 <- conf_grn_symbol[conf_grn_symbol$TF_symbol %in% tnodes, ]
# tnet2 <- conf_grn_symbol[conf_grn_symbol$Target_symbol %in% tnodes, ]
# tnet <- rbind(tnet1, tnet2)










# ###---- For later -------------------------------------------------------------------

# ## number of negative/positive regutions
# neg_num <- which(conf_grn_symbol$Direction == -1)
# pos_num <- which(conf_grn_symbol$Direction == 1)

# ## Splicing affected GRNs in cancer types
# grn_cancer <- list() 

# for(k in 1:length(all_cancer)){

#     temp <- data.table::fread(all_filesxx[k], sep='\t')
#     wh1 <- which(temp$FDR < fdr)
#     wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
#     wh <- intersect(wh1, wh2)
#     tempx <- temp[wh, ]
#     whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
#     tempy <- as.data.frame(tempx[whx,])
    
#     wh1 <- which(conf_grn_symbol$TF_symbol %in% tempy$symbol)
#     wh2 <- which(conf_grn_symbol$Target_symbol %in% tempy$symbol)

#     temp_grn <- conf_grn_symbol[union(wh1,wh2),]

#     grn_cancer[[k]] <- temp_grn

#     cat('Cancer',k,'of',length(all_cancer),'done\n')

# }


# ##--- number of negative and positive regulations ---
# neg_num_cancer <- c()
# pos_num_cancer <- c()

# for(k in 1:length(all_cancer)){

#     temp <- grn_cancer[[k]]
#     neg_num_cancer <- c(neg_num_cancer, length(which(temp$Direction == -1))/length(temp[[1]]))
#     pos_num_cancer <- c(pos_num_cancer, length(which(temp$Direction == 1))/length(temp[[1]]))

#     cat('Cancer',k,'of',length(all_cancer),'done\n')

# }


# ## mfinder downloaded from https://www.weizmann.ac.il/mcb/alon/download/network-motif-software
# ## main.c corrected --> added "extern" in the Role_members variable
# ## Once compiled, add path to ~/.profile
# save_dir <- '../data/motifs'
# if(!dir.exists(save_dir)){
#     dir.create(save_dir)
# }

# for(k in 1:length(all_cancer)){

#     temp <- grn_cancer[[k]][,c(1,2)]
#     temp$wt <- rep(1, length(temp[[1]]))
#     data.table::fwrite(temp,'../data/grn.tmp',sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
#     cmd <- paste('mfinder ../data/grn.tmp -s 3 -r 10 -f', paste0(save_dir,'/',all_cancer[k]))
#     system(cmd)

# }

