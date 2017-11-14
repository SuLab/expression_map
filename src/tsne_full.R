library(Rtsne)

load('/gpfs/group/su/lhgioia/map/data/recount/scaled/gtex_and_tcga_tpm_transpose.RData')

SRP050992.tpm <- read.table('/gpfs/group/su/lhgioia/map/data/recount/scaled/SRP050992_tpm.tsv',header=T)

row.info <- read.table('/gpfs/group/su/lhgioia/map/data/recount/row_info_genes.txt',stringsAsFactors=F,header=T)
colnames(SRP050992.tpm) <- row.info$symbol

tsne.init <- rbind(ttpm,SRP050992.tpm)

sources <- c(rep('gtex',9662),rep('tcga',11284),rep('SRP050992',459))

tsne.out <- Rtsne(as.matrix(tsne.init),perplexity=30, check_duplicates=FALSE)

tsne.out$sources <- sources

save(tsne.out,file='/gpfs/home/jbrugg/map/results/recount/tsne/tsne_tcga_gtex_SRP050992.RData')

