##############################################################################################
# Purpose: append all results in a single file
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

input_dir <- '../data/PSI_data'
fdr <- 0.05
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')

all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_files <- all_files[-4]
all_cancer <- substr(basename(all_files), 1,4)


perturbed_tfs <- '../results_new/TF_splicing/Perturbed_TF_splicing_events.xlsx'
perturbed_tfs_DBD <- '../results_new/TF_DBD/Events_perturbing_DBDs.xlsx'
perturbed_tfs_ED <- '../results_new/TF_ED/Events_perturbing_EDs.xlsx'
perturbed_tfs_survival <- '../results_new/Survival/Perturbed_TF_splicing_events_survival.xlsx'
perturbed_tfs_HD <- '../results_new/TF_HOMODIMER/Events_perturbing_HDs.xlsx'
perturbed_tfs_master <- '../results_new/MRs/Master_regulators.xlsx'

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

    openxlsx::addWorksheet(wb1, sheetName = all_cancer[k])
    openxlsx::writeData(wb1, sheet = all_cancer[k], temp1)
    openxlsx::saveWorkbook(wb1, paste0('../results_new/All_results.xlsx'), overwrite = T)

}




###---- plot for splicing events perturbing DBD, ED, both of DBD and ED, and none of DBD or ED -----------------
ed_data <- unique(data.table::fread('../data/Events_perturbing_ED.txt')$AS)
dbd_data <- unique(data.table::fread('../data/Events_perturbing_DBD.txt')$AS)
hd_data <- unique(data.table::fread('../data/Events_perturbing_HD.txt')$AS)

coding_data <- c()
for(k in 1:length(all_cancer)){
    coding_data <- union(coding_data, openxlsx::read.xlsx('../results_new/TF_coding/Events_perturbing_coding_region.xlsx',k)$AS)
}

all_splice <- c()
for(k in 1:length(all_cancer)){
    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    all_splice <- union(all_splice, tempy$as_id)
}

rem_splice <- setdiff(all_splice, union(union(union(ed_data, dbd_data), coding_data),hd_data))
no_ed_dbd_hd_splice <- setdiff(coding_data, union(union(ed_data, dbd_data),hd_data))

ed_dbd_splice <- intersect(ed_data, dbd_data)
hd_dbd_splice <- intersect(hd_data, dbd_data)
ed_hd_splice <- intersect(ed_data, hd_data)

all3_splice <- intersect(intersect(ed_data, dbd_data), hd_data)


ed_splice <- setdiff(ed_data, union(dbd_data, hd_data))
dbd_splice <- setdiff(dbd_data, union(ed_data, hd_data))
hd_splice <- setdiff(hd_data, union(ed_data,dbd_data))



pdata <- data.frame(Type=c('DNA binding domain (DBD)', 'Effector domain (ED)', 'Homodimer (HD)', 
    'Both DBD and ED', 'Both DBD and HD', 'Both ED and HD', 'DBD, ED, and HD', 'Coding but not DBD/ED/HD', 'Non-coding'), 
    countx=c(length(dbd_splice),length(ed_splice), length(hd_splice),  
        length(ed_dbd_splice),length(hd_dbd_splice),length(ed_hd_splice),length(all3_splice),
        length(no_ed_dbd_hd_splice), length(rem_splice))
    )

pdata$count <- (pdata$countx/sum(pdata$countx))*100

library(scales)
library(ggrepel)
library(tidyverse)
# blank_theme <- theme_minimal()+
#   theme(
#   axis.title.x = element_blank(),
#   axis.title.y = element_blank(),
#   panel.grid=element_blank(),
#   axis.ticks = element_blank(),
#   plot.title=element_text(size=14, face="bold")
#   )

pdata$count <- signif(pdata$count, 2)
df2 <- pdata %>% 
  mutate(csum = rev(cumsum(rev(count))), 
         pos = count/2 + lead(csum, 1),
         pos = if_else(is.na(pos), count/2, pos))

p <- ggplot(pdata, aes(x = "" , y = count, fill = fct_inorder(Type))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(count, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "% of perturbed \nsplicing events")) +
  theme_void()

# p <- ggplot(pdata, aes(x="", y=count, fill=Type))+
# geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+
# scale_fill_brewer("% of perturbed \nsplicing events") + blank_theme +
# theme(axis.text.x=element_blank())+
# geom_text_repel(aes(y = count/3 + c(0, cumsum(count)[-length(count)]), label = percent(count/100)), size=5,
#     min.segment.length = unit(0, 'lines'),)
ggsave(p,filename=paste0("../results_new/Pi_chart.png"),width=5, height=3, dpi=400)

