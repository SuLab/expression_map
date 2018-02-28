library(DESeq2)
library(BiocParallel)

register(MulticoreParam(16))

setwd('/gpfs/group/su/lhgioia/map')

allen.pca <- readRDS('results/allen/pca/allen_50_dim_4k_filtered_irlba.RDS') 

counts.mat <- as.matrix(read.table('data/allen/expected_counts_matrix.csv',sep=',',header=T))
samp.names <- rownames(counts.mat)
gene.names <- colnames(counts.mat)
counts.mat <- counts.mat * 100
counts.mat <- t(counts.mat)
counts.mat <- apply(counts.mat,2,as.integer)
rownames(counts.mat) <- gene.names
colnames(counts.mat) <- samp.names

gene.vars <- apply(counts.mat,1,var)
counts.mat <- counts.mat[gene.vars>0,]

louvain.dat <- readRDS('results/allen/clustering/louvain_pca_filtered_4k_k30.RDS')

metadata <- data.frame(cluster = louvain.dat$membership)
rownames(metadata) <- rownames(allen.pca)

metadata <- metadata[colnames(counts.mat),,drop=F]

results <- list()

for(i in unique(louvain.dat$membership)){

    cluster.vec <- rep('notGroup',ncol(counts.mat))
    cluster.vec[i == metadata$cluster] <- 'inGroup'

    pct.express.out <- round(apply(counts.mat[,cluster.vec=='notGroup',drop=F],1,function(x) sum(x > 0) / length(x)),digits=3)
    pct.express.in <- round(apply(counts.mat[,cluster.vec=='inGroup',drop=F],1,function(x) sum(x > 0) / length(x)), digits=3)

    pct.df <- data.frame(pct.express.out,pct.express.in)
    pct.max <- apply(pct.df,1,max)
    names(pct.max) <- rownames(pct.df)
    genes.use <- names(which(x=pct.max > 0.1))

    if(length(genes.use) == 0){
        print(sprintf('No genes pass min.pct threshold for cluster %d',i))
        next
    }

    metadata$assignment <- factor(cluster.vec)

    deseq.data <- DESeqDataSetFromMatrix(countData = counts.mat[genes.use,], colData = metadata, design = ~ assignment)

    deseq.data <- estimateSizeFactors(object = deseq.data)
    deseq.data <- estimateDispersions(object = deseq.data, fit = 'local')
    deseq.data <- nbinomWaldTest(object = deseq.data)
    res <- results(object = deseq.data, contrast = c('cluster','notGroup','inGroup'), alpha = 0.05, parallel=T)

    p.vals <- res$pvalue

    genes.return <- rownames(res)

    results[[i]] <- data.frame(p.vals,row.names=genes.return)
    ## return(data.frame(p.vals,row.names=genes.return))

    break
}

saveRDS(results,'results/allen/markers/louvain_deseq2_markers_cluster1.RDS')
