library('WGCNA')

enableWGCNAThreads()

setwd('/gpfs/group/su/lhgioia/map')

options(stringsAsFactors=F)

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
var.genes <- readRDS('data/allen/pca/top_5000_rotation_genes_log_point1pcthreshold.RDS')
## var.genes <- gsub('X','',var.genes)

tpm.dat <- tpm.dat[,var.genes]

mag.cnst <- as.integer(log(min(as.vector(tpm.dat)[as.vector(tpm.dat)!=0])))
dec.cnst <- exp(mag.cnst)
    
tpm.dat <- log(tpm.dat + dec.cnst) - mag.cnst

dim(tpm.dat)

powers <- c(c(1:10), seq(from = 12, to=20, by=2))

sft <- pickSoftThreshold(tpm.dat, powerVector = powers, verbose = 5,nBreaks=100)

## saveRDS(sft,'results/allen/wgcna/softThreshold.RDS')
