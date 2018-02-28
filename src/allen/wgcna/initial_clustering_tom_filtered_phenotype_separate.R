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

beta <- 2

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

TOM.all <- readRDS('data/allen/wgcna/scaled_TOM_pheno.RDS')

unmerged.labels.pheno <- list()
unmerged.colors.pheno <- list()

set.trees <- list()

for(i in 1:3){

    setTree <- hclust(as.dist(1-TOM.all[i,,]),method='average')
    set.trees[[i]] <- setTree

    unmerged.labels <- cutreeDynamic(setTree,distM=1-TOM.all[i,,],deepSplit=2,cutHeight=0.995,minClusterSize=30,pamRespectsDendro=F)
    unmerged.labels.pheno[[i]] <- unmerged.labels

    unmerged.colors <- labels2colors(unmerged.labels)
    unmerged.colors.pheno[[i]] <- unmerged.colors

}

merging.results <- list()

for(i in 1:3){

    unmerged.mes <- moduleEigengenes(expr.dat[[i]]$data,colors=unmerged.colors.pheno[[i]])

    mes <- unmerged.mes$eigengenes
    rownames(mes) <- rownames(expr.dat[[i]])
    
    me.diss <- 1-cor(mes)
    me.tree <- hclust(as.dist(me.diss),method='average')

    me.diss.thresh <- 0.25
    merge <- mergeCloseModules(expr.dat[[i]]$data,unmerged.colors.pheno[[i]],cutHeight=me.diss.thresh,verbose=3)

    merged.colors <- merge$colors
    merged.mes <- merge$newMEs
    rownames(merge$newMEs) <- rownames(expr.dat[[i]])

    color.order <- c('grey',standardColors(50))
    module.labels <- match(merged.colors,color.order)-1

    merge$labels <- module.labels

    merging.results[[i]] <- merge

}

saveRDS(merging.results,'results/allen/wgcna/merged_results_separate_phenotype.RDS')
saveRDS(set.trees,'results/allen/wgcna/gene_trees_separate_phenotype.RDS')
