## library(exactRankTests)

setwd('/gpfs/group/su/lhgioia/map/')

## allen.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
## var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')
## allen.dat <- allen.dat[,var.genes]
## saveRDS(allen.dat,'data/allen/allen_tpm_variable_genes_log_4k.RDS')

allen.pca <- readRDS('results/allen/pca/allen_50_dim_noDropouts_filtered_irlba.RDS')

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
gene.vars <- apply(tpm.dat,2,var)
tpm.dat <- tpm.dat[,gene.vars>0]

louvain.dat <- readRDS('results/allen/clustering/louvain_pca_filtered_noDropouts_k30.RDS')

cluster.vec <- louvain.dat$membership
names(cluster.vec) <- rownames(allen.pca)

cluster.vec <- cluster.vec[rownames(tpm.dat)]

## pval.mat <- matrix(0,nrow=ncol(tpm.dat),ncol=length(unique(cluster.vec)))

pval.res <- list()

for(i in unique(cluster.vec)){

    group.vec <- rep('Outside',nrow(tpm.dat))
    group.vec[cluster.vec==i] <- 'Inside'

    pct.express.out <- round(apply(tpm.dat[group.vec=='Outside',,drop=F],2,function(x) sum(x>0) / length(x)),digits=3)
    pct.express.in <- round(apply(tpm.dat[group.vec=='Inside',,drop=F],2,function(x) sum(x>0) / length(x)),digits=3) # fraction of cells expressing genes

    pct.df <- data.frame(pct.express.out,pct.express.in)
    pct.max <- apply(pct.df,1,max)
    names(pct.max) <- rownames(pct.df)
    genes.use <- names(which(x=pct.max > 0.1))

    if(length(genes.use) == 0){
        print(sprintf('No genes pass min.pct threshold for cluster %d',i))
        next
    }

    out.dat <- apply(tpm.dat[group.vec=='Outside',,drop=F],2,function(x) log(mean(expm1(x)+1)))
    in.dat <- apply(tpm.dat[group.vec=='Inside',,drop=F],2,function(x) log(mean(expm1(x)+1)))
    total.diff <- out.dat-in.dat
    genes.diff <- names(which(abs(total.diff) > 0.25))

    genes.use <- intersect(genes.use,genes.diff)
    if(length(genes.use) == 0){
        print(sprintf('No genes pass logfc.threshold for cluster %d',i))
        next
    }

    p.val <- tryCatch(
    {
        ## sapply(1:ncol(counts.mat),function(x) {return(wilcox.test(counts.mat[,x] ~ as.factor(group.vec))$p.value)})
        apply(tpm.dat[,genes.use,drop=F],2,function(x) {test = wilcox.test(x ~ as.factor(group.vec)); return(c(test$p.value,test$statistic))})
    },
        error = function(cond) {
            print(i)
            print(unique(group.vec))
            
            message(cond)
            message(sprintf('Only %d entries for cluster %d out of %d cells',sum(cluster.vec==i),i,nrow(tpm.dat)))

            break
        })
    ## cluster.markers[[i]] <- data.frame(p.val,row.names=colnames(counts.mat))
    
    names(p.val) <- genes.use

    pval.res[[i]] <- p.val

    print(sprintf('Completed cluster %d',i))

    saveRDS(pval.res,'results/allen/markers/louvain_wilcox_markers_filtered_noDropouts.RDS')
}


## saveRDS(cluster.markers,'results/allen/markers/kmeans_wilcox_markers.RDS')
## saveRDS(pval.mat,'results/allen/markers/louvain_wilcox_markers.RDS')
## 
## clusters <- clusters[rownames(allen.dat)]
##
## cluster.dat <- list()
##
## for(i in 1:length(unique(clusters))){
##
##     group.vec <- rep('Outside',nrow(allen.dat))
##     group.vec[clusters==i] <- 'Inside'
##
##     p.val <- tryCatch(
##     {
##         sapply(1:ncol(allen.dat),function(x) {return(wilcox.test(allen.dat[,x] ~ as.factor(group.vec))$p.value)})
##      
##     },
##     error = function(cond) {
##         message(cond)
##         message(sprintf('Only %d entries for cluster %d out of %d cells',sum(group.vec[clusters==i]),i,nrow(allen.dat)))
##
##         break
##     })
##  
##     cluster.dat[[i]] <- data.frame(p.val,row.names=colnames(allen.dat))
## }
##  
## saveRDS(cluster.dat,'results/allen/markers/kmeans_8_wilcox.RDS')
