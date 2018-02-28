library(RANN)
library(igraph)
library(philentropy)

setwd('/gpfs/group/su/lhgioia/map/')

## allen.pca <- readRDS('results/allen/pca/allen_50_dim_irlba.RDS')
## allen.pca <- readRDS('results/allen/pca/allen_50_dim_filtered_irlba.RDS')
allen.pca <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS') 

## allen.meta <- read.table('data/allen/allen_tsne_meta.tsv',sep='\t',header=T)
allen.meta <- read.table('data/allen/allen_tsne_meta.csv',sep=',',header=T)

louvain.results <- list()

for(reg in unique(allen.meta$sources)){

     allen.slice <- allen.pca[sapply(subset(allen.meta,sources==reg)$rnaseq_profile_id,toString),]

    neighbors <- nn2(allen.slice,k=30)
    neighbor.jacc <- 1-distance(neighbors$nn.idx,method='jaccard')    

    jacc.graph <- graph_from_adjacency_matrix(neighbor.jacc,mode="undirected",weighted=T)
    allen.louvain <- cluster_louvain(jacc.graph)

    louvain.results[[reg]] <- allen.louvain

}

saveRDS(allen.louvain,'results/allen/clustering/louvain_pca_filtered_4k_k30_sep_pheno.RDS')
## saveRDS(allen.louvain,'results/allen/clustering/louvain_50_pca.RDS')

