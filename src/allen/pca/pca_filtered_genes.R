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
## var.genes <- readRDS('results/allen/variable_genes.RDS')

var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')

allen.mat <- allen.mat[,var.genes]

mag.cnst <- as.integer(log(min(as.vector(allen.mat)[as.vector(allen.mat)!=0])))
dec.cnst <- exp(mag.cnst)
    
allen.mat <- log(allen.mat + dec.cnst) - mag.cnst

col.center <- apply(allen.mat,2,mean)

for(i in 1:ncol(allen.mat)){
    allen.mat[,i] <- allen.mat[,i] - col.center[i]
}


gene.pca <- irlba(allen.mat,nv=50,nu=0,mult=matmul)

saveRDS(gene.pca,'results/allen/pca/allen_proj_filtered_irlba_4k.RDS')

## gene.pca <- readRDS('results/recount/pca/allen_proj_all_irlba.RDS')

corr.gene.dat <- allen.mat %*% gene.pca$v
rownames(corr.gene.dat) <- rownames(allen.mat)

saveRDS(corr.gene.dat,'results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS')
