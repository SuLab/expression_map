library(recount)

setwd('/gpfs/home/jbrugg/map/data/recount/gtex/')

## for(id in list.files('.')){

## load(sprintf('%s/rse_gene.Rdata',id))
load('rse_gene_gtex.Rdata')
    ## count.dat <- count.dat[,1:10]
count.dat <- assays(rse_gene)$counts

gene.lengths <- read.table("/gpfs/home/jbrugg/map/data/recount/gene_ensembl_info.tsv", sep="\t", stringsAsFactors = F, header=T)
rownames(gene.lengths) <- gene.lengths$gene_id
rownames(count.dat) <- gene.lengths$gene_id

    ## count.dat <- cbind(count.dat, gene.lengths[,"bp_length",drop=FALSE])

scaled.dat <- count.dat/gene.lengths$bp_length

normFactors <- colSums(scaled.dat)

scaled.dat <- (10^6) * data.matrix(scaled.dat) %*% diag(1/normFactors,nrow=length(normFactors))
colnames(scaled.dat) <- colnames(count.dat)

## saveRDS(t(scaled.dat),sprintf('%s/gene_counts_tpm.RDS',id))
saveRDS(t(scaled.dat),'tpm_gtex.RDS')
## }

## setwd('/Users/Jake/Documents/Projects/Mercator/')

## count.dat <- read.table('downloads/recount/test/counts_gene.tsv',sep='\t',header=T)
## count.dat <- count.dat[,1:10]

## gene.lengths <- read.table("downloads/recount/row_info.tsv", sep="\t", stringsAsFactors = F, header=T)

## rownames(gene.lengths) <- gene.lengths$gene_id
## rownames(count.dat) <- gene.lengths$gene_id

## ## count.dat <- cbind(count.dat, gene.lengths[,"bp_length",drop=FALSE])

## scaled.dat <- count.dat/gene.lengths$bp_length

## normFactors <- colSums(scaled.dat)

## scaled.dat <- (10^6) * data.matrix(scaled.dat) %*% diag(1/normFactors)
## colnames(scaled.dat) <- colnames(count.dat)

## write.table(t(scaled.dat),'downloads/recount/scaled/tcga_tpm.tsv')

## normFactor <- NULL

## for(sample in names(count.dat)){

##     normFactor[sample] <- as.numeric(0)
##     for(gene in rownames(count.dat)){
##         normFactor[sample] <- normFactor[sample] + (count.dat[gene,sample]/gene.lengths[gene,'bp_length'])
##     }

## }


## TPM.dat <- count.dat

## for(sample in names(count.dat)){

##     TPM.dat[,sample] <- count.dat[,sample]*(10^6)/(gene.lengths[,'bp_length']*normFactor[[sample]])
## }



