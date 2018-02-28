library('WGCNA')

enableWGCNAThreads()

setwd('/gpfs/group/su/lhgioia/map')

options(stringsAsFactors=F)

allen.meta <- read.table('data/allen/allen_meta.tsv',sep='\t',header=T)

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')

tpm.dat <- tpm.dat[as.character(allen.meta$rnaseq_profile_id),var.genes]

mag.cnst <- as.integer(log(min(as.vector(tpm.dat)[as.vector(tpm.dat)!=0])))
dec.cnst <- exp(mag.cnst)
    
tpm.dat <- log(tpm.dat + dec.cnst) - mag.cnst

dim(tpm.dat)

formatted.correct <- all(rownames(tpm.dat) == allen.meta$rnaseq_profile_id)

if(!formatted.correct){
    print('Rows out of order')
    quit()
}

powers <- c(c(1:10), seq(from = 12, to=20, by=2))

for(i in unique(allen.meta$sources)){

    cluster.dat <- tpm.dat[allen.meta$sources==i,]

    print(sprintf('Soft thresholding cluster %s...',i))
    sft <- pickSoftThreshold(cluster.dat, powerVector = powers, verbose = 100)

    ## saveRDS(sft,'results/allen/wgcna/softThreshold.RDS')
}

warnings()
