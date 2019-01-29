library(Rtsne)

setwd('/gpfs/group/su/lhgioia/map')

recount.data <- readRDS('/gpfs/group/su/lhgioia/map/results/amnio_endo/pca/recount_noAmnio_50_dim_irlba_scaled.RDS')

tsne.out <- Rtsne(recount.data,perplexity=30,check_duplicates=FALSE,pca=FALSE)

rownames(tsne.out$Y) <- rownames(recount.data)

## saveRDS(tsne.out,file='/gpfs/group/su/lhgioia/map/results/recount/tsne/recount_tsne.RDS')
saveRDS(tsne.out,file='/gpfs/group/su/lhgioia/map/results/amnio_endo/tsne/recount_noAmnio_tsne_p30_scaledPCA.RDS')
