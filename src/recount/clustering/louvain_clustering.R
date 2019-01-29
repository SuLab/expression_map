library(RANN)
library(igraph)
library(philentropy)
library(sets)
library(arrangements)
library(parallel)

getSimilarity <- function(x){

    y <- nn.recount[x[1],]
    z <- nn.recount[x[2],]

    union <- length(union(z,y))
    inter <- length(intersect(z,y))

    sim <- inter/union
    if(sim < 1/15){
        sim <- 0
    }

    return(sim)
}

setwd('/gpfs/group/su/lhgioia/map/')

## allen.pca <- readRDS('results/allen/pca/allen_50_dim_irlba.RDS')
## allen.pca <- readRDS('results/allen/pca/allen_50_dim_filtered_irlba.RDS')
## allen.pca <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS') 

## recount.pca <- readRDS('results/recount/pca/recount_250_dim_noScaled_noProj_over100_noSingle.RDS')

## ## neighbors <- nn2(allen.pca,k=30)
## neighbors <- nn2(recount.pca,k=30)

## nn.recount <- read.table('results/recount/clustering/vptree_test_inds.csv',sep=',',stringsAsFactors=F,row.names=1)
nn.recount <- read.table('results/recount/clustering/recount_nn_k100.csv',sep=',',stringsAsFactors=F,row.names=1)
nn.recount <- nn.recount+1

## neighbor.sets <- list()

print('beginning')

## for(i in 1:nrow(nn.recount)){
##     neighbor.sets[[i]] <- as.set(nn.recount[i,])
## }

neighbor.jacc <- matrix(1,nrow=nrow(nn.recount),ncol=nrow(nn.recount))

## print('made sets')

no_cores <- detectCores()-1
cl <- makeCluster(no_cores)

clusterExport(cl,'nn.recount')

ind.combs <- combinations(nrow(nn.recount),2,replace=T)

weights <- parApply(cl,ind.combs,1,getSimilarity)

## i <- 1
## while(i < nrow(nn.recount)){    

##     for(j in (i+1):nrow(nn.recount)){

##         val <- set_similarity(neighbor.sets[[i]],neighbor.sets[[j]])

##         if(val < 1/15) val <- 0

##         neighbor.jacc[i,j] <- val
##         neighbor.jacc[j,i] <- val

##     }
##     i <- i+1
## }

## neighbor.jacc <- 1-distance(neighbors$nn.idx,method='jaccard')

print('calculated similarities')

stopCluster(cl)

saveRDS(weights,'results/recount/clustering/nn_sim_weights_over50_k100.RDS')

jacc.graph <- graph_from_adjacency_matrix(neighbor.jacc,mode="upper",weighted=T)
edge_attr(jacc.graph,'weight') <- weights

print('made graph')

recount.louvain <- cluster_louvain(jacc.graph)

print('finished clustering')

recount.louvain$names <- rownames(nn.recount)

## neighbor.jacc <- 1-distance(neighbors$nn.idx,method='jaccard')
## jacc.graph <- graph_from_adjacency_matrix(neighbor.jacc,mode="undirected",weighted=T)
## allen.louvain <- cluster_louvain(jacc.graph)

saveRDS(recount.louvain,'results/recount/clustering/louvain_pca_over50_noSingle_k100.RDS')
## saveRDS(allen.louvain,'results/allen/clustering/louvain_50_pca.RDS')


