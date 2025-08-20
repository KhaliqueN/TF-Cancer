##############################################################################################
# Purpose: Cancer switching events
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(pheatmap)

save_dir <- '../results_new/TF_switching'

if(dir.exists(save_dir)){
    unlink(save_dir, recursive=TRUE)
}
dir.create(save_dir, recursive=TRUE)

## TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
paired_sam <- data.table::fread('../data/cancer_paired_samples.txt')
paired_sam <- paired_sam[-4,]
##----------------------------------------------------------------------

input_dir <- '../data/PSI_data'
fdr <- 0.05
num_sams <- 0.5
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI_paired.txt', full.names=TRUE))
all_files <- all_files[-4]
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))
all_files_raw <- all_files_raw[-4]
all_cancer <- substr(basename(all_files), 1,4)

##------------------------------------------------
events_tf <- list()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    wh <- which(temp$FDR < fdr)
    tempx <- temp[wh, ]
    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- tempx[whx,]
    tempy$POSP <- tempy$POS/paired_sam$PAIRED_SAMPLES[k]
    tempy$NEGP <- tempy$NEG/paired_sam$PAIRED_SAMPLES[k]

    events_tf[[k]] <- tempy

}


###----- switching events ------------------------------------------------------------------------------------
loop1 <- length(all_cancer)-1
loop2 <- loop1+1
pos1 <- c()
neg1 <- c()
pos2 <- c()
neg2 <- c()
mdiff1 <- c()
mdiff2 <- c()
tcancer1 <- c()
tcancer2 <- c()
asid <- c()
tsymbol <- c()
tstype <- c()
for(k in 1:loop1){
    i <- k+1
    temp1 <- events_tf[[k]]
    for(j in i:loop2){
        temp2 <- events_tf[[j]]
        both_eve <- intersect(temp1$as_id, temp2$as_id)
        ## for each common event --------------------------
        for(ce in 1:length(both_eve)){
            temp1a <- temp1[temp1$as_id == both_eve[ce], ]
            temp2a <- temp2[temp2$as_id == both_eve[ce], ]
            tcancer1 <- c(tcancer1, all_cancer[k])
            tcancer2 <- c(tcancer2, all_cancer[j])
            mdiff1 <- c(mdiff1, temp1a$MEDIAN_DIFF)
            mdiff2 <- c(mdiff2, temp2a$MEDIAN_DIFF)
            asid <- c(asid, both_eve[ce])
            tsymbol <- c(tsymbol, temp1a$symbol)
            tstype <- c(tstype, temp1a$splice_type)

            pos1 <- c(pos1, temp1a$POSP)
            neg1 <- c(neg1, temp1a$NEGP)
            pos2 <- c(pos2, temp2a$POSP)
            neg2 <- c(neg2, temp2a$NEGP)
        }
    }
}

pdata <- data.frame(SYMBOL=tsymbol, AS_ID=asid, SPLICE_TYPE=tstype, CANCER1=tcancer1, CANCER2=tcancer2, MEDIAN_DIFF1=mdiff1, MEDIAN_DIFF2=mdiff2, POS1=pos1, NEG2=neg2, POS2=pos2, NEG1=neg1)


##---- Plot the number of switching events with different cutoffs ----
perc <- seq(0.5,1,0.1)
clen <- c()
for(k in 1:length(perc)){
    wha <- which(pdata$POS1 >= perc[k] & pdata$NEG2 >= perc[k])
    # pdata1 <- pdata[wh,]
    whb <- which(pdata$NEG1 >= perc[k] & pdata$POS2 >= perc[k])
    # pdata2 <- pdata[wh,]
    wh <- union(wha, whb)
    pdata1 <- pdata[wh, ]
    clen <- c(clen, length(pdata1[[1]]))
}

temp_data <- data.frame(pert=perc, clen=clen)
temp_data$pert <- as.factor(temp_data$pert)
basesize <- 10
ppx <- ggplot(data = temp_data, aes(x=pert, y=clen)) + 
geom_bar(stat="identity")+
scale_x_discrete()+
scale_y_continuous(, limits=c(0,max(temp_data$clen)+350))+
xlab("Fraction of patients with \npaired samples")+ylab("# of cases with cancer-type \nswitching splicing events")+
geom_text(aes(label=clen), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
# guides(color=guide_legend(title="Significance",ncol=1,override.aes = list(size = 2)))
guides(color='none')
ggsave(ppx,filename=paste0(save_dir, "/Switching_events_numbers.png"),width=3.5, height=3.5, dpi=400)



##---- Then show the delta PSI distribution of a chosen cutoff for the cancer pairs with most number of switching events ----
chosen_cut <- 0.8
wha <- which(pdata$POS1 >= chosen_cut & pdata$NEG2 >= chosen_cut)
# pdata1 <- pdata[wh,]
whb <- which(pdata$NEG1 >= chosen_cut & pdata$POS2 >= chosen_cut)
# pdata2 <- pdata[wh,]
wh <- union(wha, whb)
pdata1 <- pdata[wh, ]
data.table::fwrite(pdata1, paste0(save_dir,'/Switching_events.csv'))

pdata1$ID <- paste0(pdata1$SYMBOL,'_',pdata1$AS_ID,'_',pdata1$SPLICE_TYPE, '\n(',pdata1$CANCER1,'_',pdata1$CANCER2,')')
all_ids <- unique(pdata1$ID) ## 57 out of the possible 91 cancer type pairs have atleast one switching event

pdata2 <- pdata1[,c(2,3)]
all_ids_count1 <- plyr::count(pdata2)
all_ids_count1 <- all_ids_count1[order(-all_ids_count1$freq),]

# pdata3 <- pdata1[pdata1$CANCER1 == 'KICH' & pdata1$CANCER2 == 'KIRP', ]
# uasid <- unique(pdata3$AS_ID)

# for(k in 1:length(uasid)){

#     temp_kirc <- data.table::fread(all_files[6])
#     temp_kirc <- temp_kirc[temp_kirc$as_id == uasid[k], ]
# }

tolabel1 <- subset(pdata1, abs(MEDIAN_DIFF1) > 0.2)$ID
tolabel2 <- subset(pdata1, abs(MEDIAN_DIFF2) > 0.2)$ID
tolabel <- intersect(tolabel1, tolabel2)
whl <- which(pdata1$ID %in% tolabel)
pdata1$PL <- ""
pdata1$PL[whl] <- tolabel

basesize <- 10
ppx <- ggplot(data = pdata1, aes(x=MEDIAN_DIFF1, y=MEDIAN_DIFF2, label=PL)) + 
geom_point(alpha=0.8, size=1)+
scale_x_continuous()+
scale_y_continuous()+
# scale_color_manual(values=c('#1b9e77', '#d95f02', '#7570b3', '#d9d9d9'))+
xlab("Median \u0394PSI of cancer type 1")+ylab("Median \u0394PSI of cancer type 2")+
geom_text_repel( family = "Poppins",
    size = 2.5,
    min.segment.length = 0.1, 
    seed = 42, 
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = 0,
    nudge_y = 0,
    color = "black") +
theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
# guides(color=guide_legend(title="Significance",ncol=1,override.aes = list(size = 2)))
guides(color='none')
ggsave(ppx,filename=paste0(save_dir, "/Switching_events.png"),width=4, height=3.5, dpi=400)



##---- Show pathways differences ??



# ##--- scatter plot -----------------------------
# wha <- which(pdata$POS1 > 0.75 & pdata$NEG2 > 0.75)
# # pdata1 <- pdata[wh,]
# whb <- which(pdata$NEG1 > 0.75 & pdata$POS2 > 0.75 )
# # pdata2 <- pdata[wh,]
# wh <- union(wha, whb)
# pdata1 <- pdata[wh, ]

# basesize <- 10
# ppx <- ggplot(data = pdata1, aes(x=MEDIAN_DIFF1, y=MEDIAN_DIFF2)) + 
# geom_point(alpha=0.8, size=0.5)+
# scale_x_continuous()+
# scale_y_continuous()+
# # scale_color_manual(values=c('#1b9e77', '#d95f02', '#7570b3', '#d9d9d9'))+
# xlab("Median of paired sample differences \nfor cancer type 1 (cancer-normal)")+ylab("Median of paired sample differences \nfor cancer type 2 (cancer-normal)")+
# # geom_text_repel( family = "Poppins",
# #     size = 2,
# #     min.segment.length = 0.1, 
# #     seed = 42, 
# #     box.padding = 0.5,
# #     max.overlaps = Inf,
# #     arrow = arrow(length = unit(0.010, "npc")),
# #     nudge_x = 0,
# #     nudge_y = 0,
# #     color = "black") +
# theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
#     axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
#     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"))+
# # guides(color=guide_legend(title="Significance",ncol=1,override.aes = list(size = 2)))
# guides(color='none')
# ggsave(ppx,filename=paste0(save_dir, "/Switching_events.png"),width=4, height=3.5, dpi=400)

