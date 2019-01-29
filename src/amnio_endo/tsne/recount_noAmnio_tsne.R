library(Rtsne)

setwd('/gpfs/group/su/lhgioia/map')

recount.data <- readRDS('/gpfs/group/su/lhgioia/map/results/recount/pca/recount_noAmnio_50_dim_irlba.RDS')

tsne.out <- Rtsne(recount.data,perplexity=150,check_duplicates=FALSE,pca=FALSE)

rownames(tsne.out$Y) <- rownames(recount.data)

## saveRDS(tsne.out,file='/gpfs/group/su/lhgioia/map/results/recount/tsne/recount_tsne.RDS')
saveRDS(tsne.out,file='/gpfs/group/su/lhgioia/map/results/recount/tsne/recount_noAmnio_tsne_p150.RDS')
