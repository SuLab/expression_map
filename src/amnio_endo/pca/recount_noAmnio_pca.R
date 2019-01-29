library(bigmemory)
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

## tpm.proj.mat <- readRDS('data/recount/tpm_projected_recount.RDS')

## tcga.rot <- readRDS('data/recount/tcga/proj_tcga.RDS')

## print('tcga')

## gtex.rot <- readRDS('data/recount/gtex/proj_gtex.RDS')

## print('gtex')

## data.mat <- rbind(tcga.rot,gtex.rot)


## for(file.id in list.files('data/recount/project_cnts')){

##     print(file.id)

##     ## proj.tpm <- t(as.matrix(readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_tpm.RDS',file.id))))
##     ## proj.rot <- t(rot.matrix %*% proj.tpm)
##     ## rownames(proj.rot) <- colnames(proj.tpm)
##     ## saveRDS(proj.rot,sprintf('data/recount/project_cnts/%s/gene_counts_proj.RDS',file.id))

##     proj.rot <- readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_proj.RDS',file.id))

##     data.mat <- rbind(data.mat,proj.rot)

## }

## saveRDS(data.mat,'data/recount/tpm_projected_recount.RDS')

## data.mat <- readRDS('data/recount/tpm_projected_recount.RDS')

## recount.meta <- read.table('data/recount/pheno/all_recount_metasra_summarized.tsv',sep='\t',header=T,stringsAsFactors=F)

## amnio.meta <- subset(recount.meta,proj_id == 'SRP015439')
## ## amnio.mat <- data.mat[rownames(data.mat) %in% amnio.meta[,4],]
## ## saveRDS(amnio.mat,'data/amnio_endo/amnio_tpm_samples.RDS')


## #remove NAs
## ## recount.meta <- subset(recount.meta,!all(is.na(c(organ,tissue,efo_disease,doid_disease,disease_staging,efo_cl,cl_cl))))

## ## recount.meta <- subset(recount.meta,!is.na(organ) | !is.na(tissue) | !is.na(efo_disease) | !is.na(doid_disease) | !is.na(disease_staging) | !is.na(efo_cl) | !is.na(cl_cl))
## #remove amnio experiment samples
## ## recount.meta <- subset(recount.meta,proj_id != 'SRP015439')

## ## amnio.mat <- data.mat[!(rownames(data.mat) %in% recount.meta[,4]),]
## ## saveRDS(amnio.mat,'data/amnio_endo/amnio_tpm_samples.RDS')

## data.mat <- data.mat[rownames(data.mat) %in% recount.meta[,4],]

## print(nrow(data.mat))

## if(nrow(data.mat) == 0){
##     quit()
## }

## ## col.centers <- apply(data.mat,2,mean)
## ## saveRDS(col.centers,'data/recount/col_centers_noAmnio.RDS')

## col.vars <- apply(data.mat,2,var)
## saveRDS(col.vars,'data/recount/col_variances_noAmnio.RDS')

## ## below this un-comment later

## ## big.data <- as.big.matrix(data.mat,backingfile='noAmnio_recount_projTPM.bin',descriptorfile='noAmnio_recount_projTPM.desc',backingpath='/gpfs/group/su/lhgioia/map/results/recount/pca/tmp')

big.data <- attach.big.matrix('/gpfs/group/su/lhgioia/map/results/recount/pca/tmp/noAmnio_recount_projTPM.desc')

## big.data <- attach.big.matrix('/gpfs/group/su/lhgioia/map/results/recount/pca/tmp/desc')

for(i in 1:ncol(big.data)){
    big.data[,i] <- big.data[,i] - mean(big.data[,i])
    big.data[,i] <- big.data[,i] / sqrt(var(big.data[,i]))
}

print('beginning pca')

gene.pca <- irlba(big.data,nv=50,nu=0,mult=matmul)

saveRDS(gene.pca,'/gpfs/group/su/lhgioia/map/results/recount/pca/recount_proj_noAmnio_big_irlba_scaled.RDS')

corr.gene.dat <- as.matrix(big.data %*% gene.pca$v)
rownames(corr.gene.dat) <- rownames(big.data)

saveRDS(corr.gene.dat,'/gpfs/group/su/lhgioia/map/results/recount/pca/recount_noAmnio_50_dim_irlba_scaled.RDS')

