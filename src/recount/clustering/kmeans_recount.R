setwd('/gpfs/group/su/lhgioia/map')

## allen.pca <- readRDS('results/allen/pca/allen_50_dim_irlba.RDS')
## allen.pca <- readRDS('results/allen/pca/allen_50_dim_filtered_irlba.RDS')
## allen.tpm <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)

print('start')

recount.pca <- readRDS('results/recount/pca/recount_250_dim_noScaled_noProj_over50_noSingle.RDS')

print(dim(recount.pca))

kmeans.ass <- list()

for(i in 1:100){

    print(i)

    recount.kmeans <- kmeans(recount.pca,centers=i,nstart=20)

    kmeans.ass[[i]] <- recount.kmeans$cluster

    print('clustered')

    saveRDS(kmeans.ass,'results/recount/clustering/kmeans_1to100_recount_250_dim_noScaled_noProj_over50_noSingle.RDS')

    print('saved')
}

saveRDS(kmeans.ass,'results/recount/clustering/kmeans_1to100_recount_250_dim_noScaled_noProj_over50_noSingle.RDS')

