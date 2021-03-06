library(irlba)
library(bigmemory)
library(bigalgebra)
library(Rcpp)

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

setwd('/gpfs/group/su/lhgioia/map')

tcga.tpm <- as.matrix(readRDS('data/recount/tcga/tpm_tcga.RDS'))
gtex.tpm <- as.matrix(readRDS('data/recount/gtex/tpm_gtex.RDS'))

## tpm.mat <- rbind(tcga.tpm,gtex.tpm)

tpm.mat <- matrix(0,nrow=70603,ncol=58037)

cnt <- 1

tpm.mat[cnt:nrow(tcga.tpm),] <- tcga.tpm
cnt <- cnt + nrow(tcga.tpm)

tpm.rownames <- rownames(tcga.tpm)

tpm.mat[cnt:(cnt + nrow(gtex.tpm)-1),] <- gtex.tpm
cnt <- cnt + nrow(gtex.tpm)

tpm.rownames <- c(tpm.rownames,rownames(gtex.tpm))



for(file.id in list.files('data/recount/project_cnts')){

    file.tpm <- as.matrix(readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_tpm.RDS',file.id)))

    tpm.mat[cnt:(cnt + nrow(file.tpm) - 1),] <- file.tpm
    cnt <- cnt + nrow(file.tpm)

    tpm.rownames <- c(tpm.rownames,rownames(file.tpm))
    ## tpm.mat <- rbind(tpm.mat,file.tpm)

}

print(cnt)
print(dim(tpm.mat))

rownames(tpm.mat) <- tpm.rownames

## print('making big matrix')

## big.data <- as.big.matrix(tpm.mat,backingfile='all_tpm_noProj.bin',descriptorfile='all_tpm_noProj.desc',backingpath='/gpfs/group/su/lhgioia/map/results/recount/pca/tmp')

## for(i in 1:ncol(big.data)){
##     big.data[,i] - big.data[,i] - mean(big.data[,i])
## }

for(i in 1:ncol(tpm.mat)){
    tpm.mat[,i] - tpm.mat[,i] - mean(tpm.mat[,i])
}


## print('attaching big matrix')

## big.data <- attach.big.matrix('/gpfs/group/su/lhgioia/map/results/recount/pca/tmp/all_tpm_noProj.desc')

print ('beginning pca')

gene.pca <- irlba(tpm.mat,nv=250,nu=0,mult=matmul)

saveRDS(gene.pca,'/gpfs/group/su/lhgioia/map/results/recount/pca/recount_all_big_irlba_250.RDS')

## gene.pca <- readRDS('/gpfs/group/su/lhgioia/map/results/recount/pca/recount_all_big_irlba.RDS')

rot.gene.dat <- as.matrix(tpm.mat %*% gene.pca$v)
rownames(rot.gene.dat) <- rownames(tpm.mat)

saveRDS(rot.gene.dat,'/gpfs/group/su/lhgioia/map/results/recount/pca/recount_250_dim_noScaled_noProj.RDS')
