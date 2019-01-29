library(ggplot2)
library(parallel)

getDistance <- function(x,gene.vec) {

    eu.dist <- sum((x - gene.vec)^2)
    ## corr.dist <- cor(x,gene.vec)

    return(eu.dist)

}

setwd('/gpfs/group/su/lhgioia/map/')

tsne.recount <- readRDS('results/recount/tsne/recount_noAmnio_tsne.RDS')

tcga.tpm <- as.matrix(readRDS('data/recount/tcga/tpm_tcga.RDS'))
gtex.tpm <- as.matrix(readRDS('data/recount/gtex/tpm_gtex.RDS'))

data.mat <- rbind(tcga.tpm,gtex.tpm)

rm(gtex.tpm)
rm(tcga.tpm)
gc()

for(file.id in list.files('data/recount/project_cnts')){

    print(file.id)

    proj.tpm <- as.matrix(readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_tpm.RDS',file.id)))
    ## proj.rot <- t(rot.matrix %*% proj.tpm)
    ## rownames(proj.rot) <- colnames(proj.tpm)
    ## saveRDS(proj.rot,sprintf('data/recount/project_cnts/%s/gene_counts_proj.RDS',file.id))

    ## proj.rot <- readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_proj.RDS',file.id))

    data.mat <- rbind(data.mat,proj.tpm)

}

saveRDS(data.mat,'data/recount/recount_tpm_mat.RDS')

## amnio.samples <- readRDS('data/amnio_endo/amnio_tpm_samples.RDS')

recount.meta <- read.table('data/recount/pheno/all_recount_metasra_summarized.tsv',sep='\t',header=T,stringsAsFactors=F)

amnio.meta <- subset(recount.meta,proj_id == 'SRP015439')
amnio.samples <- data.mat[rownames(data.mat) %in% amnio.meta[,4],]

data.mat <- data.mat[!(rownames(data.mat) %in% amnio.meta[,4]),]

print(dim(amnio.samples))

print(dim(data.mat))

dist.mat <- matrix(0,nrow=nrow(amnio.samples),ncol=nrow(data.mat))

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

for(i in 1:nrow(amnio.samples)){

    row <- amnio.samples[i,]

    ## center.vec <- row - map.centers
    ## center.vec <- center.vec / map.scales

    ## rot.vec <- t(center.vec) %*% rotation.dat$v

    dist.vec <- parApply(cl,data.mat,1,getDistance,gene.vec=as.vector(row))
    ## dist.vec <- apply(map.dat,1,getDistance,gene.vec=as.vector(rot.vec))

    dist.mat[i,] <- dist.vec

}


colnames(dist.mat) <- rownames(map.dat)
rownames(dist.mat) <- rownames(amnio.samples)

saveRDS(dist.mat,'results/amnio_endo/euclidian_distances.RDS')


## amnio.meta <- read.table('data/recount/metadata/SRP015439.tsv',sep='\t',header=T,stringsAsFactors=F)





