setwd('/gpfs/group/su/lhgioia/map/')

allen.pca <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS')

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
gene.vars <- apply(tpm.dat,2,var)
tpm.dat <- tpm.dat[,gene.vars>0]

louvain.dat <- readRDS('results/allen/clustering/louvain_pca_filtered_4k_k30.RDS')

cluster.vec <- louvain.dat$membership
names(cluster.vec) <- rownames(allen.pca)

cluster.vec <- cluster.vec[rownames(tpm.dat)]
        
gene.dist <- list()
gene.dist[['mean']] <- list()
gene.dist[['median']] <- list()
gene.dist[['thirdQ']] <- list()
        
for(i in 1:length(unique(cluster.vec))){

    clust.dat <- tpm.dat[cluster.vec==i,,drop=F]

    clust.mean <- apply(clust.dat,2,function(x) log(mean(x+1)))
    clust.median <- apply(clust.dat,2,function(x) log(median(x+1)))
    clust.thirdQ <- apply(clust.dat,2,function(x) log(quantile(x+1,0.75)))

    gene.dist[['mean']][[i]] <- clust.mean
    gene.dist[['median']][[i]] <- clust.median
    gene.dist[['thirdQ']][[i]] <- clust.thirdQ

    print(sprintf('Completed cluster %d',i))

    saveRDS(gene.dist,'results/allen/markers/louvain_gene_distributions.RDS')

}
    
    
## saveRDS(colnames(tpm.dat),'data/allen/tmp/variable_genes_allen.RDS')

## fc.matrix <- matrix(0,nrow=length(unique(cluster.vec)),ncol=length(unique(cluster.vec)))

## fc.gene.list <- apply(tpm.dat,2,function(x) return(fc.matrix))

## fc.res <- list()

## for(i in 1:length(unique(cluster.vec))){
##     for(j in (i+1):length(unique(cluster.vec))){

##         if(i == j) break

##         first.clust.dat <- tpm.dat[cluster.vec==i,,drop=F]
##         second.clust.dat <- tpm.dat[cluster.vec==j,,drop=F]

##         ## first.clust.dat <- apply(first.clust.dat,2,function(x) log(median(expm1(x)+1)))
##         ## second.clust.dat <- apply(second.clust.dat,2,function(x) log(median(expm1(x)+1)))

##         first.clust.dat <- apply(first.clust.dat,2,function(x) log(median(x+1)))
##         second.clust.dat <- apply(second.clust.dat,2,function(x) log(median(x+1)))

##         fc.vec <- first.clust.dat - second.clust.dat
##         names(fc.vec) <- colnames(first.clust.dat)

##         fc.res[[sprintf('%d %d',i,j)]] <- fc.vec
##         fc.res[[sprintf('%d %d',j,i)]] <- -1 * fc.vec

##         print(sprintf('Completed cluster %d %d',i,j))

##         saveRDS(fc.res,'results/allen/markers/louvain_median_fold_change_pairwise.RDS')
##     }
## }
