library(Rtsne)

setwd('/gpfs/group/su/lhgioia/map/')

rot.matrix <- as.matrix(read.table('data/recount/random_proj_rotation.csv',sep=','))

print('rotation matrix read')

tcga.tpm <- t(as.matrix(readRDS('data/recount/tcga/tpm_tcga.RDS')))
tcga.rot <- t(rot.matrix %*% tcga.tpm)
rownames(tcga.rot) <- colnames(tcga.tpm)

print('tcga rotated')

gtex.tpm <- t(as.matrix(readRDS('data/recount/gtex/tpm_gtex.RDS')))
gtex.rot <- t(rot.matrix %*% gtex.tpm)
rownames(gtex.rot) <- colnames(gtex.tpm)

print('gtex rotated')

data.mat <- rbind(tcga.rot,gtex.rot)

rm(rot.matrix)
rm(gtex.tpm)
rm(gtex.rot)
rm(tcga.tpm)
rm(tcga.rot)
gc()

sources <- c(rep('tcga',11284),rep('gtex',9662))

tsne.out <- Rtsne(data.mat,perplexity=30, check_duplicates=FALSE)

tsne.out$sources <- sources

saveRDS(tsne.out,file='/gpfs/home/jbrugg/map/results/recount/tsne/tsne_tcga_gtex.RDS')
