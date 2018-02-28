library('WGCNA')

enableWGCNAThreads()

setwd('/gpfs/group/su/lhgioia/map')

options(stringsAsFactors=F)

pca.dat <- readRDS('results/allen/pca/allen_50_dim_filtered_irlba.RDS')

louvain.clusters <- readRDS('results/allen/clustering/louvain_pca_filtered_4k_k30.RDS')

cluster.ass <- louvain.clusters$membership
names(cluster.ass) <- rownames(pca.dat)

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')

tpm.dat <- tpm.dat[,var.genes]
cluster.ass <- cluster.ass[rownames(tpm.dat)]

mag.cnst <- as.integer(log(min(as.vector(tpm.dat)[as.vector(tpm.dat)!=0])))
dec.cnst <- exp(mag.cnst)
    
tpm.dat <- log(tpm.dat + dec.cnst) - mag.cnst

dim(tpm.dat)

powers <- c(c(1:10), seq(from = 12, to=20, by=2))

for(i in unique(cluster.ass)){


    cluster.dat <- tpm.dat[cluster.ass==i,]

    print(sprintf('Soft thresholding cluster %d...',i))
    sft <- pickSoftThreshold(cluster.dat, powerVector = powers, verbose = 5,nBreaks=100)

    ## saveRDS(sft,'results/allen/wgcna/softThreshold.RDS')
}

warnings()


