library(WGCNA)

enableWGCNAThreads()

setwd('/gpfs/group/su/lhgioia/map')

options(stringsAsFactors=F)

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
var.genes <- readRDS('results/allen/variable_genes.RDS')

tpm.dat <- tpm.dat[,var.genes]

mag.cnst <- as.integer(log(min(as.vector(tpm.dat)[as.vector(tpm.dat)!=0])))
dec.cnst <- exp(mag.cnst)
    
tpm.dat <- log(tpm.dat + dec.cnst) - mag.cnst

dynamicMods <- readRDS('results/allen/wgcna/dynamic_cut_tree_filtered.RDS')
dynamicColors <- labels2colors(dynamicMods)

## MEList <- moduleEigengenes(tpm.dat,dynamicColors)

## MEs <- MEList$eigengenes
## MEDiss <- 1 - cor(MEs)

## METree <- hclust(as.dist(MEDiss),method='average')

MEDissThresh <- 0.25

merge <- mergeCloseModules(tpm.dat,dynamicColors,MEDissThresh)
## mergedColors <- merge$colors
## mergedMEs <- merge$newMEs

saveRDS(merge,'results/allen/wgcna/merged_modules_filtered.RDS')

