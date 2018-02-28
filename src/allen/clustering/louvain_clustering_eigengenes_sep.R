library(RANN)
library(igraph)
library(philentropy)
library(sets)

setwd('/gpfs/group/su/lhgioia/map/')

wgcna.eigs <- readRDS('results/allen/wgcna/projected_eigengenes_separate_pheno_filtered.RDS')

eig.mat <- NULL

for(set in names(wgcna.eigs)){

    if(is.null(eig.mat)){
        eig.mat <- wgcna.eigs[[set]]$eigengenes
    }
    else{
        set.mat <- wgcna.eigs[[set]]$eigengenes
        print('set.mat before order:')
        dim(set.mat)

        set.mat <- set.mat[rownames(eig.mat),]
        print('set.mat after order:')
        dim(set.mat)

        eig.mat <- cbind(eig.mat,set.mat)
    }

    print('eig.mat dim:')
    dim(eig.mat)

}

neighbors <- nn2(eig.mat,k=30)

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

saveRDS(allen.louvain,'results/allen/clustering/louvain_sep_eigs_k30.RDS')

