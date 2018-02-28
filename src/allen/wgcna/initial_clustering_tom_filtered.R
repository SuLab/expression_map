library(WGCNA)

enableWGCNAThreads()

setwd('/gpfs/group/su/lhgioia/map')

options(stringsAsFactors=F)

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
## var.genes <- readRDS('results/allen/variable_genes.RDS')
var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')

tpm.dat <- tpm.dat[,var.genes]

mag.cnst <- as.integer(log(min(as.vector(tpm.dat)[as.vector(tpm.dat)!=0])))
dec.cnst <- exp(mag.cnst)
    
tpm.dat <- log(tpm.dat + dec.cnst) - mag.cnst

## use threshold of 6 based off dataset size

beta <- 2
adjacency.mat <- adjacency(tpm.dat,power=beta)

rm(tpm.dat)
gc()

TOM <- TOMsimilarity(adjacency.mat)
TOM <- 1 - TOM

rm(adjacency.mat)
gc()

geneTree <- hclust(as.dist(TOM), method="average")
saveRDS(geneTree,'results/allen/wgcna/intial_gene_tree_filtered.RDS')

## geneTree <- readRDS('results/allen/wgcna/intial_gene_tree_log_4k.RDS')

minModuleSize <- 30
# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = TOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
saveRDS(dynamicMods,'results/allen/wgcna/dynamic_cut_tree_log_4k.RDS')
