library(Rtsne)

setwd('/gpfs/group/su/lhgioia/map')

allen.data <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS')

nrow(allen.data)

tsne.out <- Rtsne(allen.data,perplexity=30,check_duplicates=FALSE,pca=FALSE)

rownames(tsne.out$Y) <- rownames(allen.data)

saveRDS(tsne.out,file='results/allen/tsne/allen_tsne_filtered_log_4k.RDS')
