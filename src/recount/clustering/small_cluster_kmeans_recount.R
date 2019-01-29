setwd('/gpfs/group/su/lhgioia/map')

print('start')

recount.pca <- readRDS('results/recount/pca/recount_250_dim_noScaled_noProj_over50_noSingle.RDS')

print(dim(recount.pca))

recount.kmeans <- kmeans(recount.pca,centers=30,nstart=20)

cluster.ass <- recount.kmeans$cluster

print('clustered')

cluster.counts <- table(cluster.ass)
large.clusters <- names(cluster.counts[cluster.counts > nrow(recount.pca)*0.1])

while(length(large.clusters) > 0){

    for(cluster in large.clusters){

        cluster.slice <- recount.pca[cluster.ass == cluster,]

        sub.kmeans <- kmeans(cluster.slice,centers=5,nstart=20)
        sub.cluster.ass <- sub.kmeans$cluster
        
        sub.cluster.labels <- sapply(sub.cluster.ass,function(x) sprintf('%s-%d',toString(cluster),x))

        cluster.ass[cluster.ass == cluster] <- sub.cluster.labels

    }
    
    cluster.counts <- table(cluster.ass)
    large.clusters <- names(cluster.counts[cluster.counts > nrow(recount.pca)*0.1])

}

names(cluster.ass) <- rownames(recount.pca)

saveRDS(cluster.ass,'results/recount/clustering/kmeans_recount_over50_subclustered.RDS')

print('saved')


## saveRDS(kmeans.ass,'results/recount/clustering/kmeans_1to100_recount_250_dim_noScaled_noProj_over50_noSingle.RDS')
