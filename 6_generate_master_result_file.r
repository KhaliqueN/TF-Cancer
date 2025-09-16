##############################################################################################
# Purpose: append all results in a single file
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

input_dir <- '../data/PSI_data'
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')

all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_files <- all_files[-4]
all_cancer <- substr(basename(all_files), 1,4)


perturbed_tfs <- '../results_new/TF_splicing/Perturbed_TF_splicing_events.xlsx'
perturbed_tfs_DBD <- '../results_new/TF_DBD/PTSEs_perturbing_DBDs.xlsx'
perturbed_tfs_ED <- '../results_new/TF_ED/PTSEs_perturbing_EDs.xlsx'
perturbed_tfs_FP <- '../results_new/Footprinting/Footprinting_correlation.txt'
perturbed_tfs_DP <- '../results_new/Dependency/PTSE_dependency.xlsx'

dbd_tf <- c()
ed_tf <- c()
dp_tf <- c()

for(k in 1:length(all_cancer)){
    tempx <- data.table::fread(all_files[k])
    temp1 <- openxlsx::read.xlsx(perturbed_tfs_DBD, k)
    dbd_tf <- union(dbd_tf, tempx[tempx$as_id %in% temp1$AS_ID,]$symbol)

    temp1 <- openxlsx::read.xlsx(perturbed_tfs_ED, k)
    ed_tf <- union(ed_tf, tempx[tempx$as_id %in% temp1$AS_ID,]$symbol)
}

for(k in 1:12){
    temp1 <- openxlsx::read.xlsx(perturbed_tfs_DP, k)
    temp1 <- temp1[temp1$MIN_CHRONOS < -1, ]
    if(nrow(temp1) != 0){
        dp_tf <- union(dp_tf, temp1$TF)
    }
}


temp1 <- data.table::fread(perturbed_tfs_FP)
fp_tf <- unique(temp1[temp1$FDR < 0.05, ]$TF)


direct <- union(dbd_tf, ed_tf)

indirect <- union(dp_tf, fp_tf)






dbd <- data.table::fread(perturbed_tfs_DBD)

# perturbed_tfs_survival <- '../results_new/Survival/Perturbed_TF_splicing_events_survival.xlsx'
# perturbed_tfs_HD <- '../results_new/TF_HOMODIMER/Events_perturbing_HDs.xlsx'
# perturbed_tfs_master <- '../results_new/MRs/Master_regulators.xlsx'
# perturbed_tfs_atac <- '../results_new/Footprinting/Scatter_footprint_depth_flank.csv'

wb1 <- openxlsx::createWorkbook(paste0('../results_new/All_results.xlsx'))

for(k in 1:length(all_cancer)){

    if(k < 4){
        temp1 <- openxlsx::read.xlsx(perturbed_tfs, k)
    }else{
        temp1 <- openxlsx::read.xlsx(perturbed_tfs, k+1)
    }

    temp2 <- openxlsx::read.xlsx(perturbed_tfs_DBD, k)
    tdbd <- rep(NA, length(temp1[[1]]))
    for(j in 1:length(temp2[[1]])){
        wh <- which(temp1$as_id == temp2$AS[j])
        tdbd[wh] <- temp2$DBD[j]
    }
    temp1$DBD <- tdbd

    temp2 <- openxlsx::read.xlsx(perturbed_tfs_ED, k)
    tdbd <- rep(NA, length(temp1[[1]]))
    for(j in 1:length(temp2[[1]])){
        wh <- which(temp1$as_id == temp2$AS[j])
        tdbd[wh] <- temp2$DBD[j]
    }
    temp1$ED <- tdbd

    temp2 <- openxlsx::read.xlsx(perturbed_tfs_HD, k)
    tdbd <- rep(NA, length(temp1[[1]]))
    for(j in 1:length(temp2[[1]])){
        wh <- which(temp1$as_id == temp2$AS[j])
        tdbd[wh] <- temp2$DBD[j]
    }
    temp1$HD <- tdbd

    temp2 <- openxlsx::read.xlsx(perturbed_tfs_survival, k)
    tdbd <- rep(NA, length(temp1[[1]]))
    srty <- rep(NA, length(temp1[[1]]))
    for(j in 1:length(temp2[[1]])){
        wh <- which(temp1$as_id == temp2$SE[j])
        if(!is.na(tdbd[wh])){
            tdbd[wh] <- paste(tdbd[wh],temp2$Hazard_Ratio[j],collapse=';')
            srty[wh] <- paste(srty[wh],temp2$CI_Type[j],collapse=';')
        }else{
            tdbd[wh] <- temp2$Hazard_Ratio[j]
            srty[wh] <- temp2$CI_Type[j]
        }

    }
    temp1$Survival_Type <- srty
    temp1$Survival_HR <- tdbd

    temp2 <- openxlsx::read.xlsx(perturbed_tfs_master, k)
    tdbd <- rep(NA, length(temp1[[1]]))
    for(j in 1:length(temp2[[1]])){
        wh <- which(temp1$as_id == temp2$asid[j])
        tdbd[wh] <- 'Yes'
    }
    temp1$Master_regulator <- tdbd


    temp2 <- read.csv(perturbed_tfs_atac)
    temp2 <- temp2[temp2$CANCER == all_cancer[k], ]
    tdbd <- rep(NA, length(temp1[[1]]))

    if(nrow(temp2) != 0){
        for(j in 1:length(temp2[[1]])){
            wh <- which(temp1$as_id == temp2$ASID[j])
            tdbd[wh] <- temp2$SIG[wh]
        }
    }

    temp1$ATAC_correlation <- tdbd

    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], temp1)
    openxlsx::saveWorkbook(wb1, paste0('../results_new/All_results.xlsx'), overwrite = T)

}







