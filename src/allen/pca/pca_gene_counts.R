library(irlba)

matmul <- function(A,B,transpose=FALSE){

    if(is.null(dim(B))) B <- cbind(B)

    if(is.null(dim(A))) A <- rbind(A)

    if(transpose)
        return(cbind((t(B) %*% A)[]))

    ## print(typeof(A))
    ## print(typeof(B))

    ## print(dim(A))
    ## print(dim(B))

    cbind((A %*% B)[])

}


setwd('/gpfs/group/su/lhgioia/map/')

allen.mat <- as.matrix(read.table('data/allen/tpm_matrix.csv',sep=',',header=T))

mag.cnst <- as.integer(log(min(as.vector(allen.mat)[as.vector(allen.mat)!=0])))
dec.cnst <- exp(mag.cnst)
    
allen.mat <- log(allen.mat + dec.cnst) - mag.cnst

col.center <- apply(allen.mat,2,mean)
col.var <- apply(allen.mat,2,var)

## saveRDS(col.var,'results/allen/pca/gene_variances_log.RDS')

## zero.vars <- (col.var!=0)

## saveRDS(zero.vars,'results/allen/pca/allen_variable_columns.RDS')

for(i in 1:ncol(allen.mat)){
    allen.mat[,i] <- allen.mat[,i] - col.center[i]
    ## col.center <- c(col.center,mean(big.data[,i]))
}

## allen.mat <- allen.mat[,zero.vars]

## ## seed.vec <- runif(n=ncol(allen.mat),min=-10,max=10)

gene.pca <- irlba(allen.mat,nv=500,nu=0,mult=matmul)
rownames(gene.pca$v) <- colnames(allen.mat)

saveRDS(gene.pca,'results/allen/pca/allen_all_irlba_nv500.RDS')


## ## gene.pca <- readRDS('results/recount/pca/allen_proj_all_irlba.RDS')

corr.gene.dat <- allen.mat %*% gene.pca$v
rownames(corr.gene.dat) <- rownames(allen.mat)

saveRDS(corr.gene.dat,'results/allen/pca/allen_all_nv500_irlba_projected.RDS')
