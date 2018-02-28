library(Seurat)

setwd('/gpfs/group/su/lhgioia/map/')

allen.tpm <- read.table('data/allen/tpm_matrix.csv',header=T,sep=',')

allen.seurat <- CreateSeuratObject(raw.data=allen.tpm,min.cells=3,min.genes=200)


