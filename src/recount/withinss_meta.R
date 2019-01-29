library(jsonlite)

setwd('/gpfs/group/su/lhgioia/map')

tpm.mat <- matrix(0,nrow=70603,ncol=58037)

cnt <- 1

tpm.mat[cnt:nrow(tcga.tpm),] <- tcga.tpm
cnt <- cnt + nrow(tcga.tpm)

tpm.mat[cnt:(cnt + nrow(gtex.tpm)-1),] <- gtex.tpm
cnt <- cnt + nrow(gtex.tpm)

for(file.id in list.files('data/recount/project_cnts')){

    file.tpm <- as.matrix(readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_tpm.RDS',file.id)))

    tpm.mat[cnt:(cnt + nrow(file.tpm) - 1),] <- file.tpm
    cnt <- cnt + nrow(file.tpm)
    ## tpm.mat <- rbind(tpm.mat,file.tpm)

}


print('making big matrix')

## big.data <- as.big.matrix(tpm.mat,backingfile='all_tpm_noProj.bin',descriptorfile='all_tpm_noProj.desc',backingpath='/gpfs/group/su/lhgioia/map/results/recount/pca/tmp')

for(i in 1:ncol(tpm.mat)){
    tpm.mat[,i] - tpm.mat[,i] - mean(tpm.mat[,i])
}


## recount.meta <- read.table('data/recount/metadata/all_recount_metadata_ordered.tsv',sep='\t',header=T)

#### ## read in metadata here!

meta.json <- fromJSON(metadata)

for(name in names(meta.json)){

    runs <- meta.json[[name]]
    runs <- runs[runs %in% rownames(tsne.dat)]

    tsne.dat[runs,'metalabel'] <- name

}

for(label in unique(tsne.dat$metalabel)){

    

}

cast


## for proj.id in unique(recount.meta$proj_id){

##     meta.slice <- subset(big.data,

##     tpm.mat
## }
