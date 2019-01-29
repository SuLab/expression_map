library(RANN)
library(igraph)
library(philentropy)
library(sets)

setwd('/gpfs/group/su/lhgioia/map/')

recount.pca <- readRDS('results/amnio_endo/pca/recount_noAmnio_50_dim_irlba_scaled.RDS')


kmeans.ass <- list()

for(i in 1:30){

    allen.kmeans <- kmeans(recount.pca,centers=i,nstart=20)
    kmeans.ass[[i]] <- allen.kmeans

}

saveRDS(kmeans.ass,'results/amnio_endo/clustering/kmeans_1to30_allen_filtered_4k_irlba.RDS')


