library(Rtsne)

setwd('/gpfs/group/su/lhgioia/map')

allen.data <- readRDS('results/allen/pca/allen_50_dim_filtered_irlba.RDS')

tsne.out <- Rtsne(allen.data,perplexity=30,check_duplicates=FALSE,pca=FALSE)

rownames(tsne.out$Y) <- rownames(allen.data)

saveRDS(tsne.out,file='results/allen/tsne/allen_tsne_filtered.RDS')
