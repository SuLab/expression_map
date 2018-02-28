library(WGCNA)

enableWGCNAThreads()

setwd('/gpfs/group/su/lhgioia/map')
options(stringsAsFactors=F)

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
var.genes <- readRDS('results/allen/variable_genes_log.RDS')

gsg = goodSamplesGenes(tpm.dat, verbose = 3)
tpm.dat <- tpm.dat[gsg$goodSamples,gsg$goodGenes]

if (sum(!gsg$goodGenes)>0){
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
}
if (sum(!gsg$goodSamples)>0){
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
}

mag.cnst <- as.integer(log(min(as.vector(tpm.dat)[as.vector(tpm.dat)!=0])))
dec.cnst <- exp(mag.cnst)

tpm.dat <- log(tpm.dat + dec.cnst) - mag.cnst

sampleTree <- hclust(dist(tpm.dat),method="average")

