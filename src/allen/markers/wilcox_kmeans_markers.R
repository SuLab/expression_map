setwd('/gpfs/group/su/lhgioia/map/')

## allen.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
## var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')
## allen.dat <- allen.dat[,var.genes]
## saveRDS(allen.dat,'data/allen/allen_tpm_variable_genes_log_4k.RDS')

allen.pca <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS')

counts.mat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)

louvain.dat <- readRDS('results/allen/clustering/louvain_pca_filtered_4k_k30.RDS')

cluster.vec <- louvain.dat$membership
names(cluster.vec) <- rownames(allen.pca)

cluster.vec <- cluster.vec[colnames(counts.mat)]

pval.mat <- matrix(0,nrow=ncol(counts.mat),ncol=length(unique(cluster.vec)))

for(i in unique(cluster.vec)){

    group.vec <- rep('Outside',nrow(counts.mat))
    group.vec[cluster.vec==i] <- 'Inside'

    p.val <- tryCatch(
    {
        ## sapply(1:ncol(counts.mat),function(x) {return(wilcox.test(counts.mat[,x] ~ as.factor(group.vec))$p.value)})
        apply(counts.mat,2,function(x) wilcox.test(x ~ as.factor(group.vec)$p.value))
    },
        error = function(cond) {
            print(j)
            print(unique(group.vec))
            
            message(cond)
            message(sprintf('Only %d entries for cluster %d out of %d cells',sum(cluster.vec==j),j,nrow(counts.mat)))

            break
        })
    ## cluster.markers[[i]] <- data.frame(p.val,row.names=colnames(counts.mat))

    pval.mat[,i] <- p.val
}



## saveRDS(cluster.markers,'results/allen/markers/kmeans_wilcox_markers.RDS')
saveRDS(pval.mat,'results/allen/markers/louvain_wilcox_markers.RDS')


## clusters <- clusters[rownames(allen.dat)]

## cluster.dat <- list()

## for(i in 1:length(unique(clusters))){

##     group.vec <- rep('Outside',nrow(allen.dat))
##     group.vec[clusters==i] <- 'Inside'

##     p.val <- tryCatch(
##     {
##         sapply(1:ncol(allen.dat),function(x) {return(wilcox.test(allen.dat[,x] ~ as.factor(group.vec))$p.value)})
        
##     },
##     error = function(cond) {
##         message(cond)
##         message(sprintf('Only %d entries for cluster %d out of %d cells',sum(group.vec[clusters==i]),i,nrow(allen.dat)))

##         break
##     })
    
##     cluster.dat[[i]] <- data.frame(p.val,row.names=colnames(allen.dat))
## }
    
## saveRDS(cluster.dat,'results/allen/markers/kmeans_8_wilcox.RDS')

