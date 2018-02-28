setwd('/gpfs/group/su/lhgioia/map/')

getDistance <- function(x,gene.vec) {

    eu.dist <- sum((x - gene.vec)^2)

    return(eu.dist)

}

baldwin.tpm <- read.table('data/allen/baldwin_tpm_matrix.csv',sep=',',header=T)
var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')
baldwin.tpm <- baldwin.tpm <- bandwin.tpm[,var.genes]

allen.pca <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS')

allen.params <- readRDS('pca_scaling_params.RDS')

baldwin.tpm <- log(baldwin.tpm + map.params$dec.cnst) - map.params$mag.cnst
for(gene in var.genes){

    baldwin.tpm[,gene] <- baldwin.tpm[,gene] - map.params$col.center[gene]

}



results <- list()

allen.rot <- readRDS('results/allen/pca/allen_proj_filtered_irlba_4k.RDS')


for(i in 1:nrow(baldwin.tpm)){

    row.entry <- baldwin.tpm[i,]

    row.name <- rownames(baldwin.tpm)[i]

    rot.vec <- as.matrix(row.entry) %*% allen.rot$v

    dist.vec <- apply(allen.pca,1,getDistance,gene.vec=rot.vec)

    results[[i]] <- dist.vec

}

result.frame <- data.frame(matrix(unlist(results),nrow=nrow(baldwin.tpm),byrow=T))

rownames(result.frame) <- rownames(baldwin.tpm)

saveRDS(result.frame,'results/allen/allen_distances_filtered_4k.RDS')
