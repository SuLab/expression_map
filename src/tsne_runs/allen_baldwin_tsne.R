library(Rtsne)

setwd('/gpfs/group/su/lhgioia/map')

pca.dat <- readRDS('results/allen/pca/allen_baldwin_50_dim.RDS')

sources <- pca.dat$sources

pca.dat <- subset(pca.dat,select=-c(sources))

tsne.out <- Rtsne(pca.dat,perplexity=30,check_duplicates=FALSE,pca=FALSE)

rownames(tsne.out$Y) <- rownames(pca.dat)

tsne.out$sources <- sources

saveRDS(tsne.out,file='results/allen/tsne/allen_baldwin_tsne.RDS')
