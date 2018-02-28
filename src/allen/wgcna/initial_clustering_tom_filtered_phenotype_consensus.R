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

TOM.consensus <- pmin(TOM.all[1,,],TOM.all[2,,],TOM.all[3,,])

cons.tree <- hclust(as.dist(1-TOM.consensus),method='average')

min.module.size <- 30

unmerged.labels <- cutreeDynamic(dendro=cons.tree,distM = 1 - TOM.consensus,
                                 deepSplit = 2, cutHeight=0.995,
                                 minClusterSize=min.module.size,
                                 pamRespectsDendro=F)

unmerged.colors <- labels2colors(unmerged.labels)

unmerged.mes <- multiSetMEs(expr.dat,colors=NULL,universalColors <- unmerged.colors)

cons.me.diss.unmerged <- consensusMEDissimilarity(unmerged.mes)

cons.me.tree.unmerged <- hclust(as.dist(cons.me.diss.unmerged),method='average')

merge <- mergeCloseModules(expr.dat,unmerged.labels,cutHeight=0.25,verbose=3)

merged.labels <- merge$colors
merged.colors <- labels2colors(merged.labels)

cons.mes <- merge$newMEs

results <- list()

results$cons.tree <- cons.tree
results$unmerged.labels <- unmerged.labels
results$unmerged.colors <- unmerged.colors
results$unmerged.mes <- unmerged.mes
results$cons.me.tree.unmerged <- cons.me.tree.unmerged
results$merged.labels <- merged.labels
results$merged.colors <- merged.colors
results$cons.mes <- cons.mes

saveRDS(results,'results/allen/wgcna/consensus_phenotype_clustering.RDS')
