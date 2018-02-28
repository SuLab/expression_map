setwd('/gpfs/group/su/lhgioia/map')

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)

gene.means <- apply(tpm.dat,2,mean)
gene.sd <- apply(tpm.dat,2,sd)

mag.cnst <- as.integer(log(min(as.vector(tpm.dat)[as.vector(tpm.dat)!=0])))
dec.cnst <- exp(mag.cnst)
    
tpm.log <- log(tpm.dat + dec.cnst) - mag.cnst

log.means <- apply(tpm.log,2,mean)
log.sd <- apply(tpm.log,2,sd)

var.dat <- data.frame(gene.means,log.means,gene.sd,log.sd)
rownames(var.dat) <- colnames(tpm.dat)

saveRDS(var.dat,'data/allen/tmp/gene_variances.RDS')
