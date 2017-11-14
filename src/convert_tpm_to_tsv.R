## load('/gpfs/group/su/lhgioia/map/data/recount/scaled/gtex_and_tcga_tpm_transpose.RData')

setwd('/gpfs/group/su/lhgioia/map/')

tcga.tpm <- readRDS('data/recount/tcga/tpm_tcga.RDS')
write.table(tcga.tpm,'data/recount/tcga/tpm_tcga.tsv',sep='\t')

gtex.tpm <- readRDS('data/recount/gtex/tpm_gtex.RDS')
write.table(gtex.tpm,'data/recount/tcga/tpm_gtex.tsv',sep='\t')

## tsne.init <- rbind(tcga.tpm,gtex.tpm)

print('Starting to read projects')

for(subj in list.files('data/recount/project_cnts/')){
    subj.tpm <- readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_tpm.RDS',subj))

    ## print(subj)

    write.table(subj.tpm,sprintf('data/recount/project_cnts/%s/gene_counts_tpm.tsv',sep='\t'))

    ## tsne.init <- rbind(tsne.init,subj.tpm)
}

## SRP050992.tpm <- read.table('/gpfs/group/su/lhgioia/map/data/recount/scaled/SRP050992_tpm.tsv',header=T)

## row.info <- read.table('/gpfs/group/su/lhgioia/map/data/recount/row_info_genes.txt',stringsAsFactors=F,header=T)
## colnames(SRP050992.tpm) <- row.info$symbol

## tsne.init <- rbind(ttpm,SRP050992.tpm)

print('beginning pca')

## gene.pca <- prcomp(tsne.init,retx=TRUE)

## saveRDS(gene.pca,'/gpfs/group/su/lhgioia/map/results/recount/pca/recount_all_tpm.RDS')
