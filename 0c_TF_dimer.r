##############################################################################################
# Purpose: create map of EEIs of TF homodimers
##############################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)
library(GenomicDataCommons)
library(biomaRt)
library(seqinr)
library(Rcpp)
library(XML)
library(stringr)

rproc1 <- function(temp2x, resolution, pdbid, temp_uniprot){

    temp3 <- strsplit(trimws(temp2x[5]),'[=]')[[1]]
    tempe <- trimws(temp2x[3])
    ## choose temp3 entries with at least two length because that means there is chain information
    temp3_count <- length(temp3) # only taking the proteins with continuous information

    if(temp3_count == 2){

        temp4 <- temp3[1]
        chain <- strsplit(temp4, '[/]')[[1]]
        temp5 <- temp3[2]
        temp6 <- paste0(unlist(lapply(strsplit(temp5, '[.]'),'[[',1)),'-') 
        start <- unlist(lapply(strsplit(temp6, '[-]'),'[[',1))
        end <- unlist(lapply(strsplit(temp6, '[-]'),'[[',2))
        # replacing '' value in start and end vector
        start[start == ''] <- 0 
        end[end == ''] <- 0

        tuni <- vector(mode = "list", length = 6)

        if(resolution <= 3){
            for(ii in 1:length(chain)){
                chainx <- chain[ii]
                pdbc <- paste0(tolower(pdbid), '_', chainx)
                tuni[[1]][ii] <- temp_uniprot
                tuni[[2]][ii] <- resolution
                tuni[[3]][ii] <- pdbc
                tuni[[4]][ii] <- start
                tuni[[5]][ii] <- end
                tuni[[6]][ii] <- tempe
            }
        }
        
        return(tuni)

    }else{
        return(NULL)
    }

}

rproc2 <- function(temp2x, pdbid, temp_uniprot){

    temp3 <- strsplit(trimws(temp2x[5]),'[=]')[[1]]
    tempe <- trimws(temp2x[3])
    ## choose temp3 entries with at least two length because that means there is chain information
    temp3_count <- length(temp3)

    if(temp3_count == 2){

        temp4 <- temp3[1]
        chain <- strsplit(temp4, '[/]')[[1]] 
        temp5 <- temp3[2]
        temp6 <- paste0(unlist(lapply(strsplit(temp5, '[.]'),'[[',1)),'-') 
        start <- unlist(lapply(strsplit(temp6, '[-]'),'[[',1))
        end <- unlist(lapply(strsplit(temp6, '[-]'),'[[',2))
        # replacing '' value in start and end vector
        start[start == ''] <- 0 
        end[end == ''] <- 0

        tuni <- vector(mode = "list", length = 5)

        for(ii in 1:length(chain)){
            chainx <- chain[ii]
            pdbc <- paste0(tolower(pdbid), '_', chainx)
            tuni[[1]][ii] <- temp_uniprot
            tuni[[2]][ii] <- pdbc
            tuni[[3]][ii] <- start
            tuni[[4]][ii] <- end
            tuni[[5]][ii] <- tempe
        }
        
        return(tuni)

    }else{
        return(NULL)
    }
    
}


### Save all potential homodimers present in the PDB ------------------------------------
uniprot_HS_dir <- '../../../public_data/UniProt_database_HS'
allfiles <- list.files(uniprot_HS_dir, full.names=TRUE)
uniprotid1 <- c()
genename1 <- c()
pdbchain1 <- c()
resol1 <- c()
allstart1 <- c()
allend1 <- c()
expr1 <- c()

uniprotid2 <- c()
genename2 <- c()
pdbchain2 <- c()
allstart2 <- c()
allend2 <- c()
expr2 <- c()

pdb_entry <- 0
pdb_str <- 0
pdb_xr <- 0
pdb_nmr <- 0
expr_info <- c()
pdb_files <- c()

for(k in 1:length(allfiles)){

    temp <- readLines(allfiles[k])
    tempcc <- substr(temp, 1,2)
    
    # all uniprot ids
    wha <- which(tempcc == 'AC')
    # take the canonical uniprot identifier
    temp1 <- unlist(lapply(strsplit(temp[wha[1]], '[;]'), '[[',1))
    temp11 <- lapply(strsplit(temp1,'\\s+'), '[[',2)
    temp_uniprot <- temp11[[1]]

    # all pdb
    whp <- which(temp %like% ' PDB;')

    if(length(whp) == 0){
        next # skip if no pdb antry is present
    }else{
        pdb_entry <- pdb_entry+1
        temp_pdb <- temp[whp]
        temp2 <- strsplit(temp_pdb,'[;]')
        experiment <- trimws(unlist(lapply(temp2, '[[', 3)))
        expr_info <- union(expr_info, experiment)

        ## only keeping entries with X-ray, EM or NMR
        whx <- which(experiment %in% c('X-ray', 'EM', 'NMR'))
        if(length(whx) != 0){ pdb_str <- pdb_str+1 }
        temp_pdb <- temp_pdb[whx]
        temp2 <- strsplit(temp_pdb,'[;]')
        pdbid <- trimws(unlist(lapply(temp2, '[[', 2)))
        pdb_files <- union(pdb_files,unique(pdbid))
        resolution <- substr(trimws(unlist(lapply(temp2, '[[', 4))),1,4)
        resolution[resolution == '-'] <- 100
        resolution <- as.numeric(resolution)
        whr <- which(resolution <= 3)

        if(length(whr) != 0){

            pdb_xr <- pdb_xr+1

            temp_pdbx <- temp_pdb[whr]
            temp2 <- strsplit(temp_pdbx,'[;]')

            for(i in 1:length(temp2)){
                tempxx <- rproc1(temp2[[i]], resolution[whr[i]], pdbid[whr[i]], temp_uniprot)
                if(is.null(tempxx)){next}
                uniprotid1 <- c(uniprotid1, tempxx[[1]])
                pdbchain1 <- c(pdbchain1, tempxx[[3]])
                allstart1 <- c(allstart1, tempxx[[4]])
                allend1 <- c(allend1, tempxx[[5]])
                resol1 <- c(resol1, tempxx[[2]])
                expr1 <- c(expr1, tempxx[[6]])
            }

        }else{
            ## check if the resolution is 100
            whr <- which(resolution == 100)
            if(length(whr) != 0){

                pdb_nmr <- pdb_nmr+1

                temp_pdbx <- temp_pdb[whr]
                temp2 <- strsplit(temp_pdbx,'[;]')

                for(i in 1:length(temp2)){
                    tempxx <- rproc2(temp2[[i]], pdbid[whr[i]], temp_uniprot)
                    if(is.null(tempxx)){next}
                    uniprotid2 <- c(uniprotid2, tempxx[[1]])
                    pdbchain2 <- c(pdbchain2, tempxx[[2]])
                    allstart2 <- c(allstart2, tempxx[[3]])
                    allend2 <- c(allend2, tempxx[[4]])
                    expr2 <- c(expr2, tempxx[[5]])
                }
    
                

            }
        }

    }
    cat(k,' of ', length(allfiles), ' done\n')
}

id_map1 <- data.frame(uniprotkbac=uniprotid1, pdbchain=pdbchain1, resolution=resol1, start=allstart1, end=allend1, exp=expr1)
id_map2 <- data.frame(uniprotkbac=uniprotid2, pdbchain=pdbchain2, resolution=rep(0,length(uniprotid2)), 
    start=allstart2, end=allend2, exp=expr2)

id_map1$PDBID <- unlist(lapply(strsplit(id_map1$pdbchain, '[_]'), '[[', 1))
id_map2$PDBID <- unlist(lapply(strsplit(id_map2$pdbchain, '[_]'), '[[', 1))
id_map1$CHAIN <- unlist(lapply(strsplit(id_map1$pdbchain, '[_]'), '[[', 2))
id_map2$CHAIN <- unlist(lapply(strsplit(id_map2$pdbchain, '[_]'), '[[', 2))
id_map1$ID <- paste0(id_map1$uniprotkbac,'_',id_map1$PDBID)
id_map2$ID <- paste0(id_map2$uniprotkbac,'_',id_map2$PDBID)


##-- choose the PDBIDs that result in the highest number of complexes ---
allpdbs <- plyr::count(id_map1$ID)
allpdbsu <- allpdbs[allpdbs$freq > 1, ] ## ensures to only take homodimers
id_map11 <- id_map1[id_map1$ID %in% allpdbsu[[1]], ]


## -- choose the PDBIDs that results in the highest number of complexes ---
allpdbs <- plyr::count(id_map2$ID)
allpdbsu <- allpdbs[allpdbs$freq > 1, ]
id_map21 <- id_map2[id_map2$ID %in% allpdbsu[[1]], ]

cmx_data <- rbind(id_map11, id_map21)


##---- Check which TFs are present ---
tfs <- data.table::fread('../data/filtered_TFs_curated.txt', sep='\t')
ensembl_gene_map <- data.table::fread('../data/ensembl_name_map.txt')
tfs_uniprot <- setdiff(unique(ensembl_gene_map[ensembl_gene_map$Ensembl_gene_id %in% tfs$Ensembl_Gene_ID, ]$Uniprotswissprot),'')
## 1612 of 1639 TFs mapped to uniprot IDs

cmx_data_tf <- cmx_data[cmx_data$uniprotkbac %in% tfs_uniprot, ]
cmx_data_tf$len <- as.numeric(cmx_data_tf$end)-as.numeric(cmx_data_tf$start)+1
## 273 out of the 1612 TFs are present as potential homodimers (co-crystallized) in the PDB --

# data.table::fwrite(cmx_data_tf,'../data/TF_homodimers.txt',)

##-- Filter the data to keep only one PDB ID per protein ---
## I do this because I only want to detect an interacting exon pair between as complete proteins as possible.
alluni <- unique(cmx_data_tf$uniprotkbac)
cmx_data_tf_filt <- as.data.frame(matrix(nrow=0,ncol=0))
for(k in 1:length(alluni)){
    temp <- cmx_data_tf[cmx_data_tf$uniprotkbac == alluni[k], ]
    temp1 <- temp[temp$len == max(temp$len), ]
    tempids <- gtools::mixedsort(temp1$PDBID)[1]
    temp2 <- temp1[temp1$PDBID == tempids, ]
    cmx_data_tf_filt <- rbind(cmx_data_tf_filt, temp2)
}

data.table::fwrite(cmx_data_tf_filt,'../data/homodimer_complex.txt', sep='\t', row.names=FALSE)

##--- Detect exon pairs of the potential homodimers ---



##--- Download PDB data --------
keep <- unique(cmx_data_tf_filt$PDBID)

# check if the folder exists
store_cif <- '../../../public_data/PDB_CIF_temp'

if(dir.exists(store_cif)){
  pdbids2 <- unique(unlist(lapply(strsplit(list.files(store_cif), '[.]'), '[[', 1)))
  pdbids3 <- setdiff(keep, pdbids2)
  allpdbs <- paste(pdbids3, collapse=',')
}else{
  dir.create(store_cif)
  allpdbs <- paste(keep, collapse=',')
}

if(nchar(allpdbs) != 0){
  writeLines(allpdbs,'../data/all_CIFs.txt')
  system(paste0('./batch_download.sh -f ../data/all_CIFs.txt -o ',store_cif,' -c'))
  system(paste0('gunzip ',store_cif,'/*.cif.gz'))
}



##--- Download SIFTS ---------------------------------------------------------------------------
keep <- unique(cmx_data_tf_filt$PDBID)
store_dir <- '../../../public_data/SIFTS_temp'

if(dir.exists(store_dir)){
    allfiles <- list.files(store_dir)
    allpresent <- unlist(lapply(strsplit(allfiles, '[.]'), '[[', 1))
    todownload <- setdiff(keep, allpresent)
}else{
    dir.create(store_dir)
    todownload <- keep
}

## download SIFT data
allsifts <- substr(todownload, 2,3)

for(k in 1:length(todownload)){
    output_name <- paste0(todownload[k],'.xml.gz')
    query <- paste0('https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/',allsifts[k],'/',output_name)
    cmd1 <- paste0('wget -O ',store_dir,'/',output_name,' ',query)
    cmd2 <- paste0("gunzip --force ", store_dir,'/',output_name)
    system(cmd1)
    system(cmd2)
}



##-- map uniprot positions to pdb positios using SIFTS -------------------------------------------------------
store_dir <- '../data/uniprotPDB_map'
if(dir.exists(store_dir)){
    unlink(store_dir, recursive=TRUE)
}
dir.create(store_dir)


uniprot_pdb1 <- cmx_data_tf_filt
# nu <- plyr::count(uniprot_pdb1$uniprotkbac)
pdbids1 <- unique(uniprot_pdb1$PDBID)

for(k in 1:length(pdbids1)){ 
    ##-- get the sub dataframe for this pdb ----------------------------------------
    tempdata <- uniprot_pdb1[uniprot_pdb1$PDBID == pdbids1[k], ]
    temp <- xmlParse(paste0('../../../public_data/SIFTS_temp/',pdbids1[k],'.xml'))
    residue_set <- getNodeSet(temp, "//rs:residue[@dbSource='PDBe']", "rs")
    PDBeResNum <- c()
    uniprotId <- c()
    chainId <- c()
    uniprotResId <- c()
    uniprotResNum <- c()
    for(j in 1:length(residue_set)){
        tempr <- xmlToList(residue_set[[j]])
        wh <- which(names(tempr) == 'crossRefDb')
        temprr <- tempr[wh]
        # take names of the dbsource
        tempdbs <- unlist(unname(lapply(temprr, function(x) unname(x['dbSource']))))
        # check which one has 'UniProt'
        whu <- which(tempdbs == 'UniProt')
        if(length(whu) != 0){
            # check whether PDB resolution is mapped
            whp <- which(tempdbs == 'PDB')
            mappedNum <- unname(temprr[[whp]]['dbResNum'])

            # if some pdb residue number is mapped
            if(mappedNum != 'null'){
                # store author chain id
                chainId <- c(chainId, temprr[[whp]]['dbChainId'])
                # Store author residue number
                PDBeResNum <- c(PDBeResNum, mappedNum)
                # UniProt
                uniprotId <- c(uniprotId, temprr[[whu]]['dbAccessionId'])
                uniprotResNum <- c(uniprotResNum, temprr[[whu]]['dbResNum'])
                uniprotResId <- c(uniprotResId, temprr[[whu]]['dbResName'])
            }
        }
    }

    chainids <- tempdata$CHAIN
    unids <- tempdata$uniprotkbac

    for(j in 1:length(unids)){ # for each protein chain
        # which chain to keep
        whc <- which(chainId == chainids[j])
        uniprotId1 <- unname(uniprotId[whc])
        uniprotResNum1 <- unname(uniprotResNum[whc])
        uniprotResId1 <- unname(uniprotResId[whc])
        PDBeResNum1 <- unname(PDBeResNum[whc])
        # which uniprot to keep
        whc <- which(uniprotId1 == unids[j])
        uniprotId1 <- unname(uniprotId1[whc])
        uniprotResNum1 <- unname(uniprotResNum1[whc])
        uniprotResId1 <- unname(uniprotResId1[whc])
        PDBeResNum1 <- unname(PDBeResNum1[whc])
        # save file
        Data1 <- data.frame(UNIPROT_SEQ_NUM=uniprotResNum1, UNIPROT=uniprotResId1,PDBResNumAuthor=PDBeResNum1)
        fwrite(Data1, paste0(store_dir,'/',paste0(unids[j],'_',pdbids1[k],'_',chainids[j]),'.txt'), sep='\t', quote=FALSE, row.names=FALSE)
    }
    cat('Protein ', k, ' of ', length(pdbids1), ' done\n')
}


###---- add CIF numbers ----------------------------------------------------------------------------------------
allfiles1 <- list.files('../data/uniprotPDB_map', full.names=TRUE)

store_dir <- '../data/uniprotPDB_map_final'
if(dir.exists(store_dir)){unlink(store_dir, recursive=TRUE)}
dir.create(store_dir)

## Map CIF numbers ---------------------------------
##--------------------------------------------------------------------------------------------------------------
cif_folder <- '../../../public_data/PDB_CIF_temp'
allfiles <- list.files(cif_folder, full.names=TRUE)

for(k in 1:length(allfiles1)){

    # if((k==810) | (k==930)| (k==1320)| (k==1321)| (k==1385)| (k==1386)| (k==2914)| (k==3908)){next}# 
    ## skipped because allmap dataframe has duplicates because of "letter-code" of PDB author being missed by CIF formatting --> Database entry problem
    ## Now automatically handled later in the code --> line 91

    temp <- strsplit(basename(allfiles1[k]),'[_]')[[1]]
    temp_uni <- temp[1]
    temp_pdb <- temp[2]
    temp_chain <- strsplit(temp[3],'[.]')[[1]][1]

    wh <- which(allfiles %like% temp_pdb)
    tfile <- readLines(allfiles[wh])

    # Extract start and end positions of the coordinates entries
    wh <- which(tfile == "loop_")+1
    tfile0 <- trimws(tfile[wh])
    whh1 <- which(tfile0 == "_atom_site.group_PDB")
    start <- wh[whh1]
    if(whh1 == length(wh)){
        end <- length(tfile)
    }else{
        end <- wh[whh1+1]-1-2
    }

    # Extract the coordinates part of the PDB file
    tfile <- tfile[start:end]
    lineID <- stringr::word(tfile, 1)
    wh <- which(lineID == "ATOM")

    # Extract the field entries
    whf <- setdiff(seq(1,length(tfile)), wh)
    fields <- trimws(tfile[whf])

    tfile1 <- trimws(tfile[wh])
    tfile2 <- read.table(textConnection(tfile1), sep='', colClasses = "character")

    # author seq num
    auth_seq <- which(fields == "_atom_site.auth_seq_id")

    # cif seq num
    cif_seq <- which(fields == "_atom_site.label_seq_id")


    #filter using chain
    chainPosition <- which(fields == "_atom_site.auth_asym_id") # author-based chain ID
    chain <- tfile2[[chainPosition]]
    chain[is.na(chain)] <- "NA"
    wh <- which(chain == temp_chain)
    tfile3 <- tfile2[wh, ]
    tfile4 <- unique(tfile3[, c(cif_seq,auth_seq)])


    # map to the mapping file
    if(file.info(allfiles1[k])$size != 0){ 

        allmap <- fread(allfiles1[k])
        toadd <- rep('-', length(allmap[[1]]))

        # to catch cases where PDBResNum has a character but that is not always present in the CIF file
        tempnum <- unlist(str_extract_all(as.character(allmap$PDBResNumAuthor), "[0-9]+"))

        # to catch negative numbers
        whneg <- which(substr(as.character(allmap$PDBResNumAuthor), 1, 1) == '-')
        if(length(whneg) != 0){
            tempnum[whneg] <- paste0('-',tempnum[whneg])
        }

        whi <- intersect(tempnum, tfile4[[2]])
        wh1 <- which(tempnum %in% whi)
        wh2 <- which(tfile4[[2]] %in% whi)
        if(length(wh1)!= length(wh2)){next}## skipping cases where PDBResNum has a character but that is not always present in the CIF file
        toadd[wh1] <- tfile4[[1]][wh2]

        allmap$PDBResNumCIF <- toadd
        fwrite(allmap, paste0(store_dir,'/',basename(allfiles1[k])), sep='\t', quote=FALSE, row.names=FALSE)
    }
    cat('PDB', k, 'of', length(allfiles1), 'done\n')

}



##------- create final map -------------------------------------

# all uniprot-exon mapped files
allfiles1 <- list.files('../data/uniprot_Ensembl_Exon_map', full.names=TRUE)

# all uniprot-pdb mapped files
allfiles2 <- list.files('../data/uniprotPDB_map_final', full.names=TRUE)

# store the mappings
store_dir <- '../data/uniprot_EnsemblExonPDB_map'

if(dir.exists(store_dir)){
    unlink(store_dir, recursive=TRUE)
}
dir.create(store_dir)

# store not mapped files
nomap <- c()

# for each of the genes in the uniprot_pdb_f
for(k in 1:length(allfiles1)){

    # uniprot exon map file
    tempmapUE <- data.table::fread(allfiles1[k],header=TRUE)

    # uniprot PDB map file
    wh <- which(allfiles2 %like% strsplit(basename(allfiles1[k]),'[.]')[[1]][1])
    tempfile <- allfiles2[wh]

    if(length(tempfile) != 0){ # check whether this is present in mapped files
        for(j in 1:length(tempfile)){

            tempmapUP <- data.table::fread(tempfile[j])
            tempmapUP <- unique(tempmapUP)
            # PDB author seq number entry
            temppdb <- rep('-', length(tempmapUE[[1]]))
            tempcif <- rep('-', length(tempmapUE[[1]]))

            whi <- intersect(tempmapUE$UNIPROT_SEQ_NUM, tempmapUP$UNIPROT_SEQ_NUM)
            wh <- which(tempmapUE$UNIPROT_SEQ_NUM %in% whi)
            temppdb[wh] <- tempmapUP$PDBResNumAuthor # place the pdb seq nums
            tempcif[wh] <- tempmapUP$PDBResNumCIF # place the pdb seq nums

            temp1 <- cbind(tempmapUE[,c(3,4)],tempmapUE[,c(2)])
            temp3 <- cbind(temp1,tempmapUE[,c(1)])
            colnames(temp3) <- c('EXON', 'EXON_NUM','UNIPROT','UNIPROT_SEQ_NUM')

            #-- final map ---
            temp3$PDBResNumAuthor <- temppdb
            temp3$PDBResNumCIF <- tempcif

            fwrite(temp3, paste0(store_dir,'/',basename(tempfile[j])), row.names=FALSE, sep='\t', quote=FALSE)

        }
        
    }else{
        nomap <- c(nomap, k)
        next
    }
    
    cat('Seq ', k, 'of', length(allfiles1), 'done\n')

}

fp <- unique(unlist(lapply(strsplit(list.files('../data/uniprot_EnsemblExonPDB_map'), '[_]'),'[[',1)))

##-- 241 out of 273 as potential homodimers (co-crystallized) in the PDB are mapped--




##---- Create network --------------------------------------------------------------------------------
## PDB to unweighted network function
cppFunction("List pdb2net(CharacterVector uniqueids, NumericVector lastPositions, NumericVector xcoord, NumericVector ycoord, NumericVector zcoord, int cutoff){

    int loop1 = uniqueids.size()-1;
    int loop2 = uniqueids.size();
    int esize = loop2*loop2;
    CharacterVector p1(esize);
    CharacterVector p2(esize);
    NumericVector dis(esize);
    int start = 0;
    int counter = 0;

    for(int k=0; k<loop1; k++){

        int i = k+1;

        for(int j=i; j<loop2; j++){

            int startc = lastPositions[j-1]+1;
            int endc = lastPositions[j];
            int startr = start;
            int endr = lastPositions[k];
            double mindist = 100;

            for(int x=startr; x<=endr; x++){

                double xx = xcoord[x];
                double xy = ycoord[x];
                double xz = zcoord[x];

                for(int y=startc; y<=endc; y++){

                    double yx = xcoord[y];
                    double yy = ycoord[y];
                    double yz = zcoord[y];

                    double adist = sqrt(pow((yx-xx),2)+pow((yy-xy),2)+pow((yz-xz),2));
                    if(adist < mindist){
                        mindist = adist;
                    }

                }
            }

            if(mindist <= cutoff){
                p1[counter] = uniqueids[k];
                p2[counter] = uniqueids[j];
                dis[counter] = mindist;
                counter = counter+1;
            }

        }
        start = lastPositions[k]+1;

    }

    List L = List::create(p1,p2,dis,counter);
    return L;
  
}")

##-- all mapped data ----
allfiles <- list.files('../data/uniprot_EnsemblExonPDB_map', full.names=TRUE)
afile <- cmx_data_tf_filt
temp_pdb <- unique(afile$PDBID)

pdbDirectory <- '../../../public_data/PDB_CIF_temp'
cutoff <- 6
store_dir <- paste0('../data/networks_',cutoff)
if(dir.exists(store_dir)){
    unlink(store_dir, recursive=TRUE)
}
dir.create(store_dir)

for(k in 1:length(temp_pdb)){

    afile3 <- afile[afile$PDBID == temp_pdb[k], ]

    tfile <- readLines(paste0(pdbDirectory,"/",temp_pdb[k],".cif"))
    # Extract start and end positions of the coordinates entries
    wh <- which(tfile == "loop_")+1
    tfile0 <- trimws(tfile[wh])
    whh1 <- which(tfile0 == "_atom_site.group_PDB")
    start <- wh[whh1]
    if(whh1 == length(wh)){
        end <- length(tfile)
    }else{
        end <- wh[whh1+1]-1-2
    }

    # Extract the coordinates part of the PDB file
    tfile <- tfile[start:end]
    lineID <- word(tfile, 1)
    wh <- which(lineID == "ATOM" | lineID == "HETATM")

    # Extract the field entries
    whf <- setdiff(seq(1,length(tfile)), wh)
    fields <- trimws(tfile[whf])

    tfile1 <- trimws(tfile[wh])
    tfile2 <- read.table(textConnection(tfile1), sep='', colClasses = "character")

    # chain id is author defined and not CIF defined
    chainPosition <- which(fields == "_atom_site.auth_asym_id")#_atom_site.label_asym_id
    chain <- tfile2[[chainPosition]]
    chain[is.na(chain)] <- "NA"

    #-- for each complex from this pdbid ----
    comb <- as.data.frame(combn(afile3$CHAIN, 2))

    for(j in 1:length(comb)){

        #filter using chain
        wh <- which(chain == as.character(comb[[j]][1]) | chain == as.character(comb[[j]][2]))
        tfile3 <- tfile2[wh, ]

        #filter using model number. keep only the first model
        modelPosition <- which(fields == "_atom_site.pdbx_PDB_model_num")
        mdl <- unique(tfile3[[modelPosition]])
        wh <- which(tfile3[[modelPosition]] == mdl[1])
        tfile4 <- tfile3[wh, ]

        # extract coordinates for only the heavy atoms (S, O, C, N)
        atomPosition <- which(fields == "_atom_site.type_symbol")
        lineID2 <- tfile4[[atomPosition]]
        wh <- which(lineID2 == "C" | lineID2 == "S" | lineID2 == "O" | lineID2 == "N")
        tfile6 <- tfile4[wh, ]

        # keep only ATOM coordinates
        wh <- which(tfile6[[1]] == "ATOM")
        tfile6 <- tfile6[wh,]
        seqPosition <- which(fields == "_atom_site.label_seq_id")#CIF based
        aaPosition <- which(fields == "_atom_site.label_comp_id") # get the AA symbol column

        ### transformation into single chain ###########################
        # get the receptor and ligand information
        temp <- tfile6
        chains <- temp[[chainPosition]]
        uchain <- unique(chains)
        chain1 <- which(chains %in% comb[[j]][1])
        chain2 <- which(chains %in% comb[[j]][2])

        # original sequences
        q1 <- temp[[seqPosition]][chain1]
        q2 <- temp[[seqPosition]][chain2]
        temp[chainPosition] <- rep('Z', length(temp[[1]]))
        counter <- 1
        pointer0 <- temp[[seqPosition]][1]
        new_seq <- c(1)

        for(i in 2:length(temp[[1]])){

            pointer1 <- temp[[seqPosition]][i]

            if(pointer1 != pointer0){
                counter <- counter+1
            }

            pointer0 <- pointer1
            new_seq <- c(new_seq, counter)
        }

        temp[seqPosition] <- new_seq

        # extract chain positions of the two proteins
        temp1 <- temp[chain1, ]
        temp2 <- temp[chain2, ]
        p1 <- temp1[[seqPosition]]
        p2 <- temp2[[seqPosition]]
        a1 <- temp1[[aaPosition]]
        a2 <- temp2[[aaPosition]]

        #call to create networks
        seqq <- temp[[seqPosition]]
        wh <- which(fields == "_atom_site.Cartn_x")
        xcoord <- temp[[wh]]
        wh <- which(fields == "_atom_site.Cartn_y")
        ycoord <- temp[[wh]]
        wh <- which(fields == "_atom_site.Cartn_z")
        zcoord <- temp[[wh]]

        fname <- paste0(temp_pdb[k],'_',comb[[j]][1],'_',comb[[j]][2])
        uniqueids <- unique(seqq)
        lastPositions <- length(seqq)-match(unique(seqq),rev(seqq))
        xx <- pdb2net(as.character(uniqueids), as.numeric(lastPositions), as.numeric(xcoord), as.numeric(ycoord), as.numeric(zcoord), as.numeric(cutoff[mm]))


        Data <- data.frame(x=xx[[1]][1:xx[[4]]], y=xx[[2]][1:xx[[4]]])
        fwrite(Data,paste0(store_dir,'/',fname,'.txt'), quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
        
        pp1 <- data.frame(original=q1, new=p1, AA=a1)
        pp2 <- data.frame(original=q2, new=p2, AA=a2)

        fwrite(pp1,paste0(store_dir,'/',fname,'.chain1'), quote=FALSE, sep='\t', row.names=FALSE)
        fwrite(pp2,paste0(store_dir,'/',fname,'.chain2'), quote=FALSE, sep='\t', row.names=FALSE)

    }
    cat('Chain ',k, ' of ', length(temp_pdb), ' done\n')
}



###--------------- Exon Exon contact ------------------------------------------------

##-- all TFs which were one-2-2one mapped to a transcript in Ensembl --- 2460
allfs <- unlist(lapply(strsplit(list.files('../data/uniprot_Ensembl_Exon_map'),'[.]'),'[[',1))

cmx_f <- cmx_data_tf_filt[cmx_data_tf_filt$uniprotkbac %in% allfs, ]

mapDirectory <- '../data/uniprot_EnsemblExonPDB_map'
net_dir <- paste0('../data/networks_',cutoff)
temp_pdb <- unique(cmx_f$ID)

protein1 <- c()
protein2 <- c()
uniprot1 <- c()
uniprot2 <- c()
exon1_len <- c()
exon2_len <- c()
exon1_cov <- c()
exon2_cov <- c()
ex1 <- c()
ex2 <- c()

chainn1 <- c()
chainn2 <- c()
exon1n <- c()
exon2n <- c()
iaa1 <- c()
iaa2 <- c()
notp <- c()

for(k in 1:length(temp_pdb)){
    cmx_fx <- cmx_f[cmx_f$ID == temp_pdb[k], ]
    #-- for each complex from this pdbid ----
    comb <- as.data.frame(combn(cmx_fx$CHAIN, 2))
    for(jj in 1:length(comb)){
        t1 <- tolower(cmx_fx$PDBID[1])
        u1 <- cmx_fx$uniprotkbac[1]
        c1 <- comb[[jj]][1]
        c2 <- comb[[jj]][2]
        p1 <- paste0(t1,'_',comb[[jj]][1])
        p2 <- paste0(t1,'_',comb[[jj]][2])
        cmx_name <- paste0(t1,'_',c1,'_',c2)
        temp1 <- fread(list.files(mapDirectory, pattern=paste0('^',u1,'_',t1,'_',c1), full.names=TRUE), sep='\t', header=TRUE)
        temp2 <- fread(list.files(mapDirectory, pattern=paste0('^',u1,'_',t1,'_',c2), full.names=TRUE), sep='\t', header=TRUE)
        #------- only consider the resolved residues
        temp1 <- temp1[temp1$PDBResNumCIF != '-', ]
        temp2 <- temp2[temp2$PDBResNumCIF != '-', ]
        if(nrow(temp1) == 0 | nrow(temp2) == 0){ # if at least one mapping file has no resolved residues
            notp <- c(notp, k)
            next
        }
        temp1_ex <- unique(temp1$EXON)
        temp2_ex <- unique(temp2$EXON)
        loop1 <- length(temp1_ex)
        loop2 <- length(temp2_ex)
        if(file.exists(paste0(net_dir,'/',cmx_name,'.chain1'))){ # checking whether the name of network is saved as representing the first chain1 first or the second chain first
            cmx_name <- cmx_name
        }else{
            cmx_name <- paste0(t1,'_',c2,'_',c1)
        }
        chain1 <- unique(fread(paste0(net_dir,'/',cmx_name,'.chain1'), header=TRUE))
        chain2 <- unique(fread(paste0(net_dir,'/',cmx_name,'.chain2'), header=TRUE))
        cmx_net <- fread(paste0(net_dir,'/',cmx_name,'.txt'))
        for(i in 1:loop1){
            temp11 <- temp1[temp1$EXON == temp1_ex[i], ]
            est1 <- min(as.numeric(temp11$PDBResNumCIF))
            eed1 <- max(as.numeric(temp11$PDBResNumCIF))
            seq1 <- seq(est1, eed1) # seq nums defined by CIF
            seqq1 <- unique(chain1[which(chain1$original %in% seq1),][[2]]) # get seq nums defined by me
            for(j in 1:loop2){
                temp22 <- temp2[temp2$EXON == temp2_ex[j], ]
                est2 <- min(as.numeric(temp22$PDBResNumCIF))
                eed2 <- max(as.numeric(temp22$PDBResNumCIF))
                seq2 <- seq(est2, eed2) # seq nums defined by CIF
                seqq2 <- unique(chain2[which(chain2$original %in% seq2),][[2]]) # get seq nums defined by me
                if(length(seqq1) != 0 & length(seqq2) != 0){ # if both exons have at least one amino acid resolved
                    wh1 <- which((cmx_net$V1 %in% seqq1) & (cmx_net$V2 %in% seqq2))
                    wh2 <- which((cmx_net$V2 %in% seqq1) & (cmx_net$V1 %in% seqq2))
                    wh <- union(wh1, wh2)
                    if(length(wh) != 0){ # if the two exons interact
                        if((length(wh1) != 0) & (length(wh2) != 0)){
                            temp_net1 <- cmx_net[wh1, ]
                            temp_net2 <- cmx_net[wh2, ]
                            tt1 <- temp_net2[[2]]
                            temp_net2[[2]] <- temp_net2[[1]]
                            temp_net2[[1]] <- tt1
                            temp_net <- rbind(temp_net1, temp_net2)
                        }else if(length(wh1) != 0){
                            temp_net <- cmx_net[wh1, ]
                        }else if(length(wh2) != 0){
                            temp_net <- cmx_net[wh2, ]
                            tt1 <- temp_net[[2]]
                            temp_net[[2]] <- temp_net[[1]]
                            temp_net[[1]] <- tt1
                        }
                        allnodes <- union(temp_net[[1]], temp_net[[2]])
                        ## save all nodes before updating temp_net for later use
                        ## store the PDBResNumCIF info -------
                        ch1 <- c()
                        for(m in 1:length(temp_net[[1]])){
                            ch1 <- c(ch1, chain1[which(chain1$new == temp_net[[1]][m]),]$original)
                        }
                        ch2 <- c()
                        for(m in 1:length(temp_net[[2]])){
                            ch2 <- c(ch2, chain2[which(chain2$new == temp_net[[2]][m]),]$original)
                        }
                        ## ---- remove the PDBResNumCIF for which there are no UNIPROT-based AA positions
                        sch1 <- c()
                        for(ii in 1:length(ch1)){
                            xx <- ch1[ii]
                            xxx <- which(temp1$PDBResNumCIF == xx)
                            if(length(xxx) == 0){
                                sch1 <- c(sch1, xx)
                            }
                        }
                        ## ---- remove the PDBResNumCIF for which there are no UNIPROT-based AA positions
                        sch2 <- c()
                        for(ii in 1:length(ch2)){
                            xx <- ch2[ii]
                            xxx <- which(temp2$PDBResNumCIF == xx)
                            if(length(xxx) == 0){
                                sch2 <- c(sch2, xx)
                            }
                        }
                        ## convert sch from CIF to new IDs
                        ssch1 <- c()
                        for(m in 1:length(sch1)){
                            ssch1 <- c(ssch1, chain1[which(chain1$original == sch1[m]),]$new)
                        }
                        ssch2 <- c()
                        for(m in 1:length(sch2)){
                            ssch2 <- c(ssch2, chain2[which(chain2$original == sch2[m]),]$new)
                        }
                        tokeep1 <- setdiff(temp_net[[1]], ssch1)
                        tokeep2 <- setdiff(temp_net[[2]], ssch2)
                        wh1 <- which(temp_net[[1]] %in% tokeep1)
                        wh2 <- which(temp_net[[2]] %in% tokeep2)
                        wh <- intersect(wh1, wh2)
                        temp_net <- temp_net[wh, ]
                        ##---- now catch whether new temp_net is empty ----
                        if(nrow(temp_net) != 0){
                            ## store the PDBResNumCIF info
                            ch1 <- c()
                            for(m in 1:length(temp_net[[1]])){
                                ch1 <- c(ch1, chain1[which(chain1$new == temp_net[[1]][m]),]$original)
                            }
                            ch2 <- c()
                            for(m in 1:length(temp_net[[2]])){
                                ch2 <- c(ch2, chain2[which(chain2$new == temp_net[[2]][m]),]$original)
                            }
                            tt1 <- c()
                            for(ii in 1:length(ch1)){
                                xx <- ch1[ii]
                                xxx <- which(temp1$PDBResNumCIF == xx)
                                tt1 <- c(tt1, temp1[xxx,]$UNIPROT_SEQ_NUM)
                            }                    
                            tt2 <- c()
                            for(ii in 1:length(ch2)){
                                xx <- ch2[ii]
                                xxx <- which(temp2$PDBResNumCIF == xx) 
                                tt2 <- c(tt2, temp2[xxx,]$UNIPROT_SEQ_NUM)
                            }
                            ## store for actaul AA interactions
                            exon1n <- c(exon1n, rep(temp1_ex[i], length(tt1)))
                            exon2n <- c(exon2n, rep(temp2_ex[j], length(tt2)))
                            iaa1 <- c(iaa1, tt1)
                            iaa2 <- c(iaa2, tt2)
                            chainn1 <- c(chainn1, rep(p1, length(tt1)))
                            chainn2 <- c(chainn2, rep(p2, length(tt2)))
                        } 
                        # store info
                        protein1 <- c(protein1, p1)
                        protein2 <- c(protein2, p2)
                        exon1_len <- c(exon1_len, length(seqq1)) # length based on resolved 3D structure
                        exon2_len <- c(exon2_len, length(seqq2)) # length based on resolved 3D structure
                        exon1_cov <- c(exon1_cov, length(intersect(allnodes, seqq1)))
                        exon2_cov <- c(exon2_cov, length(intersect(allnodes, seqq2)))
                        ex1 <- c(ex1, temp1_ex[i])
                        ex2 <- c(ex2, temp2_ex[j])
                        uniprot1 <- c(uniprot1, u1)
                        uniprot2 <- c(uniprot2, u1)
                    }else{ # if no interaction
                        # store info
                        protein1 <- c(protein1, p1)
                        protein2 <- c(protein2, p2)
                        uniprot1 <- c(uniprot1, u1)
                        uniprot2 <- c(uniprot2, u1)
                        exon1_len <- c(exon1_len, length(seqq1)) # length based on resolved 3D structure
                        exon2_len <- c(exon2_len, length(seqq2)) # length based on resolved 3D structure
                        exon1_cov <- c(exon1_cov, 0)
                        exon2_cov <- c(exon2_cov, 0)
                        ex1 <- c(ex1, temp1_ex[i])
                        ex2 <- c(ex2, temp2_ex[j])
                    }
                }
            }
        }
    }
    cat('Complex', k, ' of ', length(temp_pdb), ' done\n')
}


all_exon_pairs <- data.frame(protein1=uniprot1,protein2=uniprot2,chain1=protein1, chain2=protein2, exon1=ex1, exon2=ex2, exon1_length=as.numeric(exon1_len),
    exon2_length=as.numeric(exon2_len), exon1_coverage=exon1_cov, exon2_coverage=exon2_cov)

all_exon_pairs$exon1_coverage_percent <- round((all_exon_pairs$exon1_coverage/all_exon_pairs$exon1_length)*100, 2)
all_exon_pairs$exon2_coverage_percent <- round((all_exon_pairs$exon2_coverage/all_exon_pairs$exon2_length)*100, 2)
all_exon_pairs$jaccard_percent <- round(((all_exon_pairs$exon1_coverage+all_exon_pairs$exon2_coverage)/(all_exon_pairs$exon1_length+all_exon_pairs$exon2_length))*100, 2)
fwrite(all_exon_pairs, paste0('../data/exon_pairs',cutoff,'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

## exact interacting AA info
aa_data <- data.frame(chain1=chainn1, chain2=chainn2, exon1=exon1n, exon2=exon2n, exon1_PDBResNumCIF=iaa1, exon2_PDBResNumCIF=iaa2)
fwrite(aa_data, paste0('../data/aa_interactions',cutoff,'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

## distribution of the edges based on the jaccard ########
temp_exon <- fread( paste0('../data/exon_pairs',cutoff,'.txt'), sep='\t', header=TRUE)

int_exon <- temp_exon[temp_exon$exon1_coverage > 0, ]
nint_exon <- temp_exon[temp_exon$exon1_coverage == 0, ]
fwrite(int_exon, paste0('../data/int_exon_pairs',cutoff,'.txt'), sep='\t', row.names=FALSE, quote=FALSE)



##--- map the interacting exon information on to the Exon mapping file ----
wh1 <- which(int_exon$exon1_coverage >= 5)
wh2 <- which(int_exon$exon2_coverage >= 5)
wh <- intersect(wh1, wh2)
int_exonx <- int_exon[wh, ]
alltfs <- unique(int_exonx$protein1)

store_dir <- '../data/uniprot_Ensembl_Exon_mapx'
if(dir.exists(store_dir)){
    unlink(store_dir, recursive=TRUE)
}
dir.create(store_dir)

##--- all uniprot-exon mapped files
allfiles1 <- list.files('../data/uniprot_Ensembl_Exon_map', full.names=TRUE)
##--- 156 (based on the EEI definition on >= 5 residue-residue contacts) out of the 241 co-crystallized TF homodimers have at least one EEI

for(k in 1:length(allfiles1)){

    temp_map <- data.table::fread(allfiles1[k])
    hmd <- rep('0', length(temp_map[[1]]))
    temp <- int_exon[int_exon$protein1 == strsplit(basename(allfiles1[k]),'[.]')[[1]][1], ]

    if(nrow(temp) != 0){
        ##-- there could be multiple EEI instances because of multiple chains in a pdb
        ig <- igraph::graph_from_data_frame(temp[,c(5,6)], directed = FALSE)
        igg <- igraph::simplify(ig, remove.loops=FALSE, remove.multiple=TRUE)
        temp_eei <- as.data.frame(igraph::as_edgelist(igg))
        affected_exons <- union(temp_eei[[1]], temp_eei[[2]])
        wh <- which(temp_map$EXON %in% affected_exons)
        hmd[wh] <- '1'
    }
    
    temp_map$HOMODIMER <- hmd
    data.table::fwrite(temp_map,paste0(store_dir,'/',basename(allfiles1[k])), sep='\t',row.names=FALSE, quote=FALSE)
    cat('TF', k,'of',length(allfiles1),'done\n')
}




