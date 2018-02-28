setwd('/gpfs/group/su/lhgioia/map')

## allen.pca <- readRDS('results/allen/pca/allen_50_dim_irlba.RDS')
## allen.pca <- readRDS('results/allen/pca/allen_50_dim_filtered_irlba.RDS')
## allen.tpm <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)

allen.pca <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS')

kmeans.ass <- list()

for(i in 1:20){

    allen.kmeans <- kmeans(allen.pca,centers=i,nstart=20)

    kmeans.ass[[i]] <- allen.kmeans$cluster

}

saveRDS(kmeans.ass,'results/allen/clustering/kmeans_1to20_allen_filtered_4k_irlba.RDS')

