setwd('/gpfs/group/su/lhgioia/map/')

getDistance <- function(x,gene.vec) {

    eu.dist <- sum((x - gene.vec)^2)

    return(eu.dist)

}

baldwin.tpm <- read.table('data/allen/baldwin_tpm_matrix.csv',sep=',',header=T)

allen.pca <- readRDS('results/allen/pca/allen_50_dim_irlba.RDS')

allen.centers <- readRDS('data/allen/allen_center.RDS')

results <- list()

allen.rot <- readRDS('results/allen/pca/allen_proj_all_irlba.RDS')

used.cols <- readRDS('results/allen/pca/allen_variable_columns.RDS')

for(i in 1:nrow(baldwin.tpm)){

    row.entry <- baldwin.tpm[i,]
    row.entry <- row.entry[,used.cols]

    row.name <- rownames(baldwin.tpm)[i]

    row.entry <- row.entry - allen.centers

    rot.vec <- as.matrix(row.entry) %*% allen.rot$v

    dist.vec <- apply(allen.pca,1,getDistance,gene.vec=rot.vec)

    results[[i]] <- dist.vec

}

result.frame <- data.frame(matrix(unlist(results),nrow=nrow(baldwin.tpm),byrow=T))

rownames(result.frame) <- rownames(baldwin.tpm)

saveRDS(result.frame,'results/allen/allen_distances.RDS')
