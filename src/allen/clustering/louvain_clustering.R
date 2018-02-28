library(RANN)
library(igraph)
library(philentropy)
library(sets)

setwd('/gpfs/group/su/lhgioia/map/')

## allen.pca <- readRDS('results/allen/pca/allen_50_dim_irlba.RDS')
## allen.pca <- readRDS('results/allen/pca/allen_50_dim_filtered_irlba.RDS')
allen.pca <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS') 

## neighbors <- nn2(allen.pca,k=30)
neighbors <- nn2(allen.pca,k=30)

neighbor.sets <- list()

for(i in 1:nrow(neighbors$nn.idx)){

    neighbor.sets[[i]] <- as.set(neighbors$nn.idx[i,])

}

neighbor.jacc <- matrix(1,nrow=nrow(neighbors$nn.idx),ncol=nrow(neighbors$nn.idx))

i <- 1
while(i < nrow(neighbors$nn.idx)){    

    for(j in (i+1):nrow(neighbors$nn.idx)){

        val <- set_similarity(neighbor.sets[[i]],neighbor.sets[[j]])

        if(val < 1/15) val <- 0

        neighbor.jacc[i,j] <- val
        neighbor.jacc[j,i] <- val

    }
    i <- i+1
}

## neighbor.jacc <- 1-distance(neighbors$nn.idx,method='jaccard')

jacc.graph <- graph_from_adjacency_matrix(neighbor.jacc,mode="undirected",weighted=T)

allen.louvain <- cluster_louvain(jacc.graph)



## neighbor.jacc <- 1-distance(neighbors$nn.idx,method='jaccard')

## jacc.graph <- graph_from_adjacency_matrix(neighbor.jacc,mode="undirected",weighted=T)

## allen.louvain <- cluster_louvain(jacc.graph)

saveRDS(allen.louvain,'results/allen/clustering/louvain_pca_filtered_4k_k30.RDS')
## saveRDS(allen.louvain,'results/allen/clustering/louvain_50_pca.RDS')

