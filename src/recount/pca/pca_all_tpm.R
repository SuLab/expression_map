library(bigmemory)
## library(bigpca)
library(bigalgebra)
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

## rot.matrix <- as.matrix(read.table('data/recount/random_proj_rotation.csv',sep=','))

## print('rotation matrix read')

tcga.tpm <- as.matrix(readRDS('data/recount/tcga/tpm_tcga.RDS'))
## tcga.rot <- t(rot.matrix %*% tcga.tpm)
## rownames(tcga.rot) <- colnames(tcga.tpm)
## saveRDS(tcga.rot,'data/recount/tcga/proj_tcga.RDS')

## tcga.rot <- readRDS('data/recount/tcga/proj_tcga.RDS')

print('tcga')

gtex.tpm <- as.matrix(readRDS('data/recount/gtex/tpm_gtex.RDS'))
## ## gtex.rot <- t(rot.matrix %*% gtex.tpm)
## ## rownames(gtex.rot) <- colnames(gtex.tpm)
## ## saveRDS(gtex.rot,'data/recount/gtex/proj_gtex.RDS')

## gtex.rot <- readRDS('data/recount/gtex/proj_gtex.RDS')

print('gtex')

data.mat <- rbind(tcga.tpm,gtex.tpm)

rm(gtex.tpm)

rm(tcga.tpm)

gc()

print('cleanup')

for(file.id in list.files('data/recount/project_cnts')){

    print(file.id)

    proj.tpm <- t(as.matrix(readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_tpm.RDS',file.id))))
    ## proj.rot <- t(rot.matrix %*% proj.tpm)
    ## rownames(proj.rot) <- colnames(proj.tpm)
    ## saveRDS(proj.rot,sprintf('data/recount/project_cnts/%s/gene_counts_proj.RDS',file.id))

    ## proj.rot <- readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_proj.RDS',file.id))

    data.mat <- rbind(data.mat,proj.tpm)

}

## ## saveRDS(data.mat,'data/recount/tpm_projected_recount.RDS')

## print('beginning pca')
## ## rm(rot.matrix)

big.data <- as.big.matrix(data.mat,backingfile='all_tpm.bin',descriptorfile='all_tpm.desc',backingpath='/gpfs/group/su/lhgioia/map/results/recount/pca/tmp/')

## big.data <- attach.big.matrix('/gpfs/group/su/lhgioia/map/results/recount/pca/tmp/all_tpm.desc')

## col.center <- c()

for(i in 1:ncol(big.data)){
    big.data[,i] <- big.data[,i] - mean(big.data[,i])
    ## col.center <- c(col.center,mean(big.data[,i]))
}

## col.center <- apply(data.mat,2,mean)
## col.var <- apply(data.mat,2,var)

## rm(data.mat)
## gc()

print('beginning pca')

## gene.pca <- irlba(big.data,nv=50,nu=0,mult=matmul,center=col.center)

gene.pca <- irlba(big.data,nv=50,nu=0,mult=matmul)

## gene.pca <- big.PCA(big.data,n.cores=20)

## gene.pca <- prcomp(data.mat,retx=TRUE)

saveRDS(gene.pca,'/gpfs/group/su/lhgioia/map/results/recount/pca/recount_proj_all_big_irlba.RDS')

## corr.gene.dat <- PC.correct(gene.pca,dir='/gpfs/group/su/lhgioia/map/results/recount/pca/tmp',num.pca=50,n.cores=20)

corr.gene.dat <- as.matrix(big.data %*% gene.pca$v)
rownames(corr.gene.dat) <- rownames(big.data)

saveRDS(corr.gene.dat,'/gpfs/group/su/lhgioia/map/results/recount/pca/recount_50_dim_irlba.RDS')
