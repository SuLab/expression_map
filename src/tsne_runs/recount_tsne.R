library(Rtsne)

setwd('/gpfs/group/su/lhgioia/map')

recount.data <- readRDS('/gpfs/group/su/lhgioia/map/results/recount/pca/recount_50_dim_irlba.RDS')

tsne.out <- Rtsne(recount.data,perplexity=30,check_duplicates=FALSE,pca=FALSE)

rownames(tsne.out$Y) <- rownames(recount.data)

saveRDS(tsne.out,file='/gpfs/group/su/lhgioia/map/results/recount/tsne/recount_tsne.RDS')
