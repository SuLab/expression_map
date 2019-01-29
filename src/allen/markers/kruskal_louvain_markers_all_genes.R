library(onewaytests)

setwd('/gpfs/group/su/lhgioia/map/')

## allen.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
## var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')
## allen.dat <- allen.dat[,var.genes]
## saveRDS(allen.dat,'data/allen/allen_tpm_variable_genes_log_4k.RDS')

allen.pca <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS')

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
gene.vars <- apply(tpm.dat,2,var)
## var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')
tpm.dat <- tpm.dat[,gene.vars>0]
## tpm.dat <- tpm.dat[,var.genes,drop=F]
## tpm.dat <- readRDS('data/allen/allen_tpm_variable_genes_log_4k.RDS')

louvain.dat <- readRDS('results/allen/clustering/louvain_pca_filtered_4k_k30.RDS')

cluster.vec <- louvain.dat$membership
names(cluster.vec) <- rownames(allen.pca)

cluster.vec <- cluster.vec[rownames(tpm.dat)]

## pval.mat <- matrix(0,nrow=ncol(tpm.dat),ncol=length(unique(cluster.vec)))

## pval.res <- list()

## first.clust.dat <- tpm.dat[cluster.vec==i,,drop=F]
## second.clust.dat <- tpm.dat[cluster.vec==j,,drop=F]

## samp.dat <- rbind(first.clust.dat,second.clust.dat)

## group.vec <- rep('Outside',nrow(samp.dat))
## group.vec[1:nrow(first.clust.dat)] <- 'Inside'

pct.express <- round(apply(tpm.dat,2,function(x) sum(x>0) / length(x)),digits=3)
## pct.express.in <- round(apply(second.clust.dat,2,function(x) sum(x>0) / length(x)),digits=3)

## pct.df <- data.frame(pct.express.out,pct.express.in)
## pct.max <- apply(pct.df,1,max)
## names(pct.max) <- rownames(pct.df)
genes.use <- names(which(x=pct.express > 0.1))

if(length(genes.use) == 0){
    print(sprintf('No genes pass min.pct threshold for cluster %d',i))
    next
}

## out.dat <- apply(first.clust.dat,2,function(x) log(mean(expm1(x)+1)))
## in.dat <- apply(second.clust.dat,2,function(x) log(mean(expm1(x)+1)))
## total.diff <- out.dat-in.dat
## genes.diff <- names(which(abs(total.diff) > 0.25)) # NOTE: this is using the not log of the threshold. aka using 1.2x threshold for the logfc.

## genes.use <- intersect(genes.use,genes.diff)
## if(length(genes.use) == 0){
##     print(sprintf('No genes pass logfc.threshold for cluster %d',i))
##     next
## }

gene.tests <- tryCatch(
{
    ## sapply(1:ncol(counts.mat),function(x) {return(wilcox.test(counts.mat[,x] ~ as.factor(group.vec))$p.value)})
    ## apply(tpm.dat[,genes.use,drop=F],2,function(x) (x,as.factor(cluster.vec)))
    apply(tpm.dat[,genes.use,drop=F],2,function(x) {dat = data.frame(x=x,group=as.factor(cluster.vec)); return(oneway.test(x~group,dat))})
},
error = function(cond) {
    ## print(i)
    print(unique(group.vec))
    
    message(cond)
    ## message(sprintf('Only %d entries for cluster %d out of %d cells',sum(cluster.vec==i),i,nrow(samp.dat)))

    break
})
## cluster.markers[[i]] <- data.frame(p.val,row.names=colnames(counts.mat))

names(gene.tests) <- genes.use

## pval.res[[sprintf('%d %d',i,j)]] <- data.frame(pvalue = sapply(gene.tests, function(x) x$p.value),statistic = sapply(gene.tests,function(x) x$statistic))
## c(gene.tests$p.value,gene.tests$statistic)

## print(sprintf('Completed cluster %d %d',i,j))

saveRDS(gene.tests,'results/allen/markers/louvain_welchanova_markers_all_genes.RDS')


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
