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

baldwin.mat <- as.matrix(read.table('data/allen/baldwin_tpm_matrix.csv',sep=',',header=T))

sources <- c(rep('allen',nrow(allen.mat)),rep('baldwin',nrow(baldwin.mat)))

data.mat <- rbind(allen.mat,baldwin.mat)

mag.cnst <- as.integer(log(min(as.vector(data.mat)[as.vector(data.mat)!=0])))
dec.cnst <- exp(mag.cnst)
    
data.mat <- log(data.mat + dec.cnst) - mag.cnst

col.center <- apply(data.mat,2,mean)
col.var <- apply(data.mat,2,var)
zero.vars <- (col.var!=0)

for(i in 1:ncol(data.mat)){
    data.mat[,i] <- data.mat[,i] - col.center[i]
    ## col.center <- c(col.center,mean(big.data[,i]))
}

## col.var <- apply(data.mat,2,var)

data.mat <- data.mat[,zero.vars]

gene.pca <- irlba(data.mat,nv=50,nu=0,mult=matmul)

saveRDS(gene.pca,'results/allen/pca/allen_baldwin_irlba.RDS')

## gene.pca <- readRDS('results/allen/pca/allen_baldwin_irlba.RDS')

corr.gene.dat <- as.data.frame(data.mat %*% gene.pca$v)

rownames(corr.gene.dat) <- rownames(data.mat)

corr.gene.dat$sources <- sources

saveRDS(corr.gene.dat,'results/allen/pca/allen_baldwin_50_dim.RDS')
