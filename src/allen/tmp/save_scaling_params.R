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

var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')

allen.filtered <- allen.mat[,var.genes]

mag.cnst <- as.integer(log(min(as.vector(allen.filtered)[as.vector(allen.filtered)!=0])))
dec.cnst <- exp(mag.cnst)
    
allen.mat <- log(allen.mat + dec.cnst) - mag.cnst

col.center <- apply(allen.mat,2,mean)

scaling.params <- list()

scaling.params$mag.cnst <- mag.cnst
scaling.params$dec.cnst <- dec.cnst
scaling.params$col.center <- col.center

saveRDS(scaling.params,'data/allen/pca_scaling_params.RDS')
