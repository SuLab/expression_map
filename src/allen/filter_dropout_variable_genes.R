library(WGCNA)

logVMR <- function(x,mag.cnst){
    return(log2(var(2^(x + mag.cnst))/mean(2^(x + mag.cnst))))
}

expMean <- function(x,mag.cnst, dec.cnst){
    return(log2(mean(2^(x + mag.cnst)-dec.cnst)+1))
}

setwd('/gpfs/group/su/lhgioia/map')

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)

mag.cnst <- as.integer(log2(min(as.vector(tpm.dat)[as.vector(tpm.dat)!=0])))
dec.cnst <- 2^(mag.cnst)
    
tpm.dat <- log2(tpm.dat + dec.cnst) - mag.cnst

gsg <- goodSamplesGenes(tpm.dat,verbose=100)
## tpm.dat <- tpm.dat[,gsg$goodGenes]

## zero.prop <- apply(tpm.dat,2,function(x) sum(x>1) / nrow(tpm.dat))

## saveRDS(zero.prop,'results/allen/genes_zero_entries.RDS')

zero.prop <- readRDS('results/allen/genes_zero_entries.RDS')

allen.pca <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS')
louvain.dat <- readRDS('results/allen/clustering/louvain_pca_filtered_4k_k30.RDS')

cluster.vec <- louvain.dat$membership
names(cluster.vec) <- rownames(allen.pca)

cluster.vec <- cluster.vec[rownames(tpm.dat)]

## zero.prop.list <- list()

## for(i in unique(cluster.vec)){
    
##     in.dat <- tpm.dat[cluster.vec==i,,drop=F]

##     zero.prop <- apply(in.dat,2,function(x) sum(x>1) / nrow(in.dat))

##     zero.prop.list[[i]] <- zero.prop

## }

## saveRDS(zero.prop.list,'results/allen/gene_zero_entries_cluster.RDS')

## gene.means <- apply(tpm.dat,2,expMean,mag.cnst=mag.cnst,dec.cnst=dec.cnst)
## gene.disp <- apply(tpm.dat,2,logVMR,mag.cnst=mag.cnst)

## gene.means[is.na(gene.means)] <- 0
## gene.disp[is.na(gene.disp)] <- 0

## names(gene.means) <- names(gene.disp) <- colnames(tpm.dat)

## new.gene.means <- gene.means[is.finite(gene.means) & is.finite(gene.disp)]
## new.gene.disp <- gene.disp[is.finite(gene.means) & is.finite(gene.disp)]

## gene.means <- new.gene.means
## gene.disp <- new.gene.disp

## data.mean.bin <- cut(gene.means,20)
## names(data.mean.bin) <- names(gene.means)

## disp.means <- tapply(gene.disp,data.mean.bin,mean)
## disp.sd <- tapply(gene.disp,data.mean.bin,sd)

## gene.disp.scaled <- (gene.disp - disp.means[as.numeric(data.mean.bin)]) / disp.sd[as.numeric(data.mean.bin)]

## names(gene.disp.scaled) <- names(data.mean.bin)

## gene.df <- data.frame(gene.means,gene.disp,gene.disp.scaled)
## rownames(gene.df) <- names(data.mean.bin)

## saveRDS(gene.df,'results/allen/binned_gene_dispersion_log_adj.RDS')

gene.df <- readRDS('results/allen/binned_gene_dispersion_log_adj.RDS')

## gene.df <- subset(gene.df,

low.cutoff <- quantile(gene.means,probs=c(0.05))

## passed.cutoff <- names(gene.means)[which((gene.disp.scaled > 0.5) & (gene.means > 0.0125))]
                       
passed.cutoff <- names(gene.means)[which((gene.disp.scaled > 0.5) & (gene.means > low.cutoff) & (zero.prop)]

length(passed.cutoff)

saveRDS(passed.cutoff,'results/allen/variable_genes_log.RDS')
