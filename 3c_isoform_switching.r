##############################################################################################
# Purpose: Detect which TF isoforms switch across patients
##############################################################################################

rm(list=ls())
library(data.table)
library(ggplot2)
library(pheatmap)

save_dir <- '../results'
tf_ensemb_map <- as.data.frame(data.table::fread('../data/TF_ensembl_uniprot.txt', sep='\t'))

## TFs -------------------
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
##----------------------------------------------------------------------

input_dir <- '../data/PSI_data'
nsam <- 10
fdr <- 0.05
diff <- 0.1
all_files <- gtools::mixedsort(list.files(input_dir, pattern='*filtered_PSI.txt', full.names=TRUE))
all_files_raw <- gtools::mixedsort(list.files(input_dir, pattern='^PSI_download', full.names=TRUE))
all_cancer <- substr(basename(all_files), 1,4)


##--- pairs of events of the same gene median difference (between cancer and normal) distribution ----
p1 <- c()
p2 <- c()
md1 <- c()
md2 <- c()
cancer <- c()

for(k in 1:length(all_cancer)){

    temp <- data.table::fread(all_files[k], sep='\t')
    ## save number of samples ----
    samples <- c(samples, length(which(colnames(temp) %like% 'TCGA')))
    wh1 <- which(temp$FDR < fdr)
    wh2 <- which(abs(temp$MEDIAN_DIFF) > diff)
    wh <- intersect(wh1, wh2)
    tempx <- temp[wh, ]

    whx <- which(toupper(tempx$symbol) %in% tfs$Gene_Symbol) ## number of AS events concerning TFs
    tempy <- as.data.frame(tempx[whx,])
    wh <- which(colnames(tempy) %like% '_Norm')
    tempy1 <- tempy[, -wh]
    wh <- which(colnames(tempy1) %like% 'TCGA')

    ugenes <- unique(tempy1$symbol)

    for(j in 1:length(ugenes)){

        tempg <- tempy1[tempy1$symbol == ugenes[j], ]
        if(nrow(tempg) > 1){
            loop1 <- nrow(tempg)-1
            loop2 <- nrow(tempg)
            for(i in 1:loop1){
                m <- i+1
                for(ii in m:loop2){
                    p1 <- c(p1, paste0(tempg[i,]$symbol[1],'_',tempg[i,]$as_id,'_',tempg[i,]$splice_type))
                    p2 <- c(p2, paste0(tempg[ii,]$symbol[1],'_',tempg[ii,]$as_id,'_',tempg[ii,]$splice_type))
                    md1 <- c(md1, tempg[i, ]$MEDIAN_DIFF)
                    md2 <- c(md2, tempg[ii, ]$MEDIAN_DIFF)
                    cancer <- c(cancer, all_cancer[k])
                }
            }
        }
    }

}

tempdata <- data.frame(CANCER=cancer, AS1=p1, AS2=p2, MD1=md1, MD2=md2)

##--- scatter plot ----
wh1 <- which(tempdata$MD1 < 0)
wh2 <- which(tempdata$MD2 > 0)
wha <- intersect(wh1, wh2)
wh1 <- which(tempdata$MD1 > 0)
wh2 <- which(tempdata$MD2 < 0)
whb <- intersect(wh1, wh2)
wh <- union(wha, whb)
tempdatax <- tempdata[wh,]

for(k in 1:length(tempdatax[[1]])){
    if(tempdatax$MD1[k] > 0){
        tx <- tempdatax$MD1[k]
        rx <- tempdatax$AS1[k]
        tempdatax$MD1[k] <- tempdatax$MD2[k]
        tempdatax$AS1[k] <- tempdatax$AS2[k]
        tempdatax$AS2[k] <- rx
        tempdatax$MD2[k] <- tx
    }
}

p <- ggplot(tempdatax, aes(MD1, MD2, color=CANCER)) + 
geom_point(alpha=0.6)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="Difference between median cancer \n and normal PSI value", limits=c(-1,0)) +
scale_y_continuous(name="Difference between median cancer \n and normal PSI value", limits=c(0,1)) +
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0(save_dir,"/Event_pairs_median_PSI_diff.png"),width=5, height=3, dpi=400)


##--- Which of the event pairs are also have affected DBDs -----
ttdata <- tempdata
store_val <- function(y1){
    if(nrow(y1) == 0){
        flag <- 0
    }else{
        y2 <- y1[y1$DBD != '-', ]
        if(nrow(y2) < 2){
            flag <- 0
        }else{
            flag <- unique(y2$DBD)
        }
    }
    return(paste(flag,collapse=','))
}


##-- contradictory event pairs ---
ttdata <- tempdata
wh1 <- which(ttdata$MD1 < 0 & ttdata$MD2 > 0)
wh2 <- which(ttdata$MD2 < 0 & ttdata$MD1 > 0)
ttdata <- ttdata[union(wh1, wh2), ]
pdata <- plyr::count(ttdata[[1]])
p <- ggplot(pdata, aes(x, freq)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of contradictory event pairs", limits=c(0,340)) +
geom_text(aes(label=freq), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Contradictory_event_pairs_all.png"),width=3.5, height=3, dpi=400)


##--- Which of the contradictory event pairs are also have affected DBDs -----
ttdata <- tempdata

odbd1 <- c()
odbd2 <- c()
tcancer <- c()
as1 <- c()
as2 <- c()
cor1 <- c()
cor2 <- c()
counter <- 0
counterx <- 0

for(k in 1:length(all_cancer)){

    tempxx <- ttdata[ttdata$CANCER == all_cancer[k],]
    wh1 <- which(tempxx$MD1 < 0)
    wh2 <- which(tempxx$MD2 > 0)
    wha <- intersect(wh1, wh2)
    wh1 <- which(tempxx$MD1 > 0)
    wh2 <- which(tempxx$MD2 < 0)
    whb <- intersect(wh1, wh2)
    wh <- union(wha, whb)
    tempyy <- tempxx[wh, ]
    counterx <- counterx+length(tempyy[[1]])

    for(j in 1:length(tempyy[[1]])){
        txx <- strsplit(tempyy[j,]$AS1,'[_]')[[1]][1]
        teves <- c(strsplit(tempyy[j,]$AS1,'[_]')[[1]][2], strsplit(tempyy[j,]$AS2,'[_]')[[1]][2])
        tst <- c(strsplit(tempyy[j,]$AS1,'[_]')[[1]][3], strsplit(tempyy[j,]$AS2,'[_]')[[1]][3])
        tuni <- tf_ensemb_map[tf_ensemb_map$HGNC_symbol == txx, ]$Uniprotswissprot[1]
        if(is.na(tuni)){
            counter <- counter+1
            next
        } ## some gene symbols of TFs were not mapped to Ensembl

        # if("RI" %in% tst){
        #     next
        # }
        tmap <- as.data.frame(data.table::fread(paste0('../data/uniprot_Ensembl_Exon_map_DBD_AS/',tuni,'_',all_cancer[k],'.txt')))
        
        for(i in 1:length(teves)){
            if(tst[i] == 'AP'){
                y1 <- tmap[tmap$AP == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }   
            }else if(tst[i] == 'AT'){
                y1 <- tmap[tmap$AT == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }  
            }else if(tst[i] == 'ES'){
                y1 <- tmap[tmap$ES == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }  
            }else if(tst[i] == 'AD'){
                y1 <- tmap[tmap$AD == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }  
            }else if(tst[i] == 'AA'){
                y1 <- tmap[tmap$AA == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }  
            }else if(tst[i] == 'ME'){
                y1 <- tmap[tmap$ME == teves[i], ]
                ff <- store_val(y1)
                if(i == 1){
                    odbd1 <- c(odbd1, ff)
                }else{
                    odbd2 <- c(odbd2, ff)
                }  
            }else{
                if(i == 1){
                    odbd1 <- c(odbd1, 0)
                }else{
                    odbd2 <- c(odbd2, 0)
                } 
            }
        }

        as1 <- c(as1, tempyy[j,]$AS1)
        as2 <- c(as2, tempyy[j,]$AS2)
        tcancer <- c(tcancer, all_cancer[k])
        cor1 <- c(cor1, tempyy[j,]$MD1)
        cor2 <- c(cor2, tempyy[j,]$MD2)
    }

}

tdata <- data.frame(CANCER=tcancer, AS1=as1, AS2=as2, DBD1=odbd1, DBD2=odbd2, MD1=cor1, MD2=cor2)
data.table::fwrite(tdata, paste0('../data/contradictory_event_pair_DBD.txt'), row.names=FALSE, quote=FALSE, sep='\t')
wh1 <- which(tdata$DBD1 != 0)
wh2 <- which(tdata$DBD2 != 0)
wh <- union(wh1, wh2)
tdatax <- tdata[wh,]
## 357 cases where at least one one event affects the DBD
wh1a <- setdiff(wh1, wh2)
wh2a <- setdiff(wh2, wh1)
wh <- union(wh1a, wh2a)
tdatay <- tdata[wh, ]
## 340 of the 357 are those cases where one event affects the DBD and the other does not

pdata <- plyr::count(tdatay[[1]])
p <- ggplot(pdata, aes(x, freq)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="Cancer type") + 
scale_y_continuous(name="# of contradictory event pairs", limits=c(0,50)) +
geom_text(aes(label=freq), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Contradictory_event_pairs_DBD.png"),width=3.5, height=3, dpi=400)


##--- Contradictory splicing events occuring in multiple cancer types --------------------------------
combs <- list() ## store all combinations
for(k in 1:length(all_cancer)){
    combs[[k]] <- combn(all_cancer, k)
}

##-- create undirected graphs --
cancer_gr <- list()
for(k in 1:length(all_cancer)){
    temp <- tdatay[tdatay$CANCER == all_cancer[k],]
    temp1 <- temp[,c(2,3)]
    cancer_gr[[k]] <- igraph::graph_from_data_frame(temp1, directed=FALSE)
}

num_events_combs <- list()
num_events_combs_counts <- c()
for(k in 1:length(all_cancer)){

    temp_comb <- as.data.frame(combs[[k]])
    loop1 <- length(temp_comb)
    loop2 <- length(temp_comb[[1]])
    temp_unn <- igraph::make_empty_graph(n = 0, directed = FALSE)

    for(i in 1:loop1){
        temp_ovl <- cancer_gr[[which(all_cancer == temp_comb[[i]][1])]]
        if(loop2 > 1){
            for(j in 2:loop2){
                temp_ovl <- igraph::intersection(temp_ovl, cancer_gr[[which(all_cancer == temp_comb[[i]][j])]], keep.all.vertices = TRUE)
            }
        }
        temp_unn <- igraph::union(temp_unn, temp_ovl)
    }

    num_events_combs[[k]] <- temp_unn
    num_events_combs_counts <- c(num_events_combs_counts, igraph::ecount(temp_unn))
    cat('Cancer',k,'of',length(all_cancer),'done\n')
}

pdata <- data.frame(cancer=as.factor(seq(1,length(all_cancer))), count=num_events_combs_counts)
p <- ggplot(pdata, aes(cancer, count)) + 
geom_bar(stat="identity")+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_discrete(name="# of cancer types") + 
scale_y_continuous(name="# of contradictory event pairs", limits=c(0,180)) +
geom_text(aes(label=count), hjust=0, vjust=0, angle=75, size=3)+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(fill=guide_legend(title="Percent Spliced In",ncol=1))#guides(fill='none')
ggsave(p,filename=paste0(save_dir,"/Contradictory_event_pairs_DBD_overlap.png"),width=7, height=3, dpi=400)

