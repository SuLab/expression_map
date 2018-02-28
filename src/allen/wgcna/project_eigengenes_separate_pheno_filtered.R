library(WGCNA)

enableWGCNAThreads()

setwd('/gpfs/group/su/lhgioia/map')

options(stringsAsFactors=F)

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')
allen.meta <- read.table('data/allen/allen_meta.tsv',sep='\t',header=T)

tpm.dat <- tpm.dat[as.character(allen.meta$rnaseq_profile_id),var.genes]

mag.cnst <- as.integer(log(min(as.vector(tpm.dat)[as.vector(tpm.dat)!=0])))
dec.cnst <- exp(mag.cnst)
    
tpm.dat <- log(tpm.dat + dec.cnst) - mag.cnst

expr.dat <- vector(mode='list',length=3)

cnt <- 1

for(i in unique(allen.meta$sources)){

    dat.slice <- tpm.dat[allen.meta$sources==i,]

    expr.dat[[cnt]] <- list(data=dat.slice)
    ## colnames(expr.dat[[cnt]]$data) <- colnames(dat.slice)
    ## rownames(expr.dat[[cnt]]$data) <- rownames(dat.slice)

    cnt <- cnt + 1
}

names(expr.dat) <- unique(allen.meta$sources)

merging.results <- readRDS('results/allen/wgcna/merged_results_separate_phenotype.RDS')

names(merging.results) <- names(expr.dat)

proj.results <- list()

for(pheno in names(expr.dat)){

    mes <- moduleEigengenes(tpm.dat,colors=merging.results[[pheno]]$colors)
    rownames(mes$eigengenes) <- rownames(tpm.dat)
    proj.results[[pheno]] <- mes

}

names(proj.results) <- names(expr.dat)

saveRDS(proj.results,'results/allen/wgcna/projected_eigengenes_separate_pheno_filtered.RDS')

