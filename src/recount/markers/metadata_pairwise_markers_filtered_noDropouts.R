library(jsonlite)

setwd('/gpfs/group/su/lhgioia/map/')

# make tpm matrix
 
tcga.tpm <- as.matrix(readRDS('data/recount/tcga/tpm_tcga.RDS'))
gtex.tpm <- as.matrix(readRDS('data/recount/gtex/tpm_gtex.RDS'))

tpm.dat <- matrix(0,nrow=70603,ncol=58037)

cnt <- 1

tpm.dat[cnt:nrow(tcga.tpm),] <- tcga.tpm
cnt <- cnt + nrow(tcga.tpm)

tpm.rownames <- rownames(tcga.tpm)

tpm.dat[cnt:(cnt + nrow(gtex.tpm)-1),] <- gtex.tpm
cnt <- cnt + nrow(gtex.tpm)

tpm.rownames <- c(tpm.rownames,rownames(gtex.tpm))

for(file.id in list.files('data/recount/project_cnts')){

    file.tpm <- as.matrix(readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_tpm.RDS',file.id)))

    tpm.dat[cnt:(cnt + nrow(file.tpm) - 1),] <- file.tpm
    cnt <- cnt + nrow(file.tpm)

    tpm.rownames <- c(tpm.rownames,rownames(file.tpm))
    ## tpm.dat <- rbind(tpm.dat,file.tpm)

}

print(cnt)
print(dim(tpm.dat))

rownames(tpm.dat) <- tpm.rownames
colnames(tpm.dat) <- colnames(tcga.tpm)
# filter non-variable genes

gene.vars <- apply(tpm.dat,2,var)
tpm.dat <- tpm.dat[,gene.vars>0]

print(dim(tpm.dat))

metadata.df <- read.table('data/recount/metadata/mesh_uberon_subjtree_table.csv',sep=',',stringsAsFactors=F,header=T)
rownames(metadata.df) <- metadata.df$id

metadata.info <- metadata.df['A',]
metadata.list <- fromJSON(metadata.info$termTree)

interesting.systems <- c('Digestive System','Nervous System','Endocrine System','Cardiovascular System','Musculoskeletal System','Respiratory System','Hemic and Immune Systems','Integumentary System')

## louvain.dat <- readRDS('results/allen/clustering/louvain_pca_filtered_noDropouts_k30.RDS')

## cluster.vec <- louvain.dat$membership
## names(cluster.vec) <- rownames(allen.pca)

## cluster.vec <- cluster.vec[rownames(tpm.dat)]

## pval.mat <- matrix(0,nrow=ncol(tpm.dat),ncol=length(unique(cluster.vec)))

pval.res <- list()

for(i in 1:length(interesting.systems)){
    for(j in (i+1):length(interesting.systems)){

        if(i == j | j > length(interesting.systems)){
            next
        }

        first.clust.samps <- intersect(rownames(tpm.dat),metadata.list[[interesting.systems[i]]])
        second.clust.samps <- intersect(rownames(tpm.dat),metadata.list[[interesting.systems[j]]])

        first.clust.dat <- tpm.dat[first.clust.samps,,drop=F]
        second.clust.dat <- tpm.dat[second.clust.samps,,drop=F]

        samp.dat <- rbind(first.clust.dat,second.clust.dat)

        group.vec <- rep('Outside',nrow(samp.dat))
        group.vec[1:nrow(first.clust.dat)] <- 'Inside'

        pct.express.out <- round(apply(first.clust.dat,2,function(x) sum(x>0) / length(x)),digits=3)
        pct.express.in <- round(apply(second.clust.dat,2,function(x) sum(x>0) / length(x)),digits=3)

        pct.df <- data.frame(pct.express.out,pct.express.in)
        pct.max <- apply(pct.df,1,max)
        names(pct.max) <- rownames(pct.df)
        genes.use <- names(which(x=pct.max > 0.1))

        if(length(genes.use) == 0){
            print(sprintf('No genes pass min.pct threshold for cluster %d',i))
            next
        }

        ## out.dat <- apply(first.clust.dat,2,function(x) log(mean(expm1(x)+1)))
        ## in.dat <- apply(second.clust.dat,2,function(x) log(mean(expm1(x)+1)))

        out.dat <- apply(first.clust.dat,2,function(x) log2(mean(x+1)))
        in.dat <- apply(second.clust.dat,2,function(x) log2(mean(x+1)))
        
        total.diff <- out.dat-in.dat
        genes.diff <- names(which(abs(total.diff) > 1)) # NOTE: this is using the not log of the threshold. aka using 1.2x threshold for the logfc.

        genes.use <- intersect(genes.use,genes.diff)
        if(length(genes.use) == 0){
            print(sprintf('No genes pass logfc.threshold for cluster %d',i))
            next
        }

        gene.tests <- tryCatch(
        {
            ## sapply(1:ncol(counts.mat),function(x) {return(wilcox.test(counts.mat[,x] ~ as.factor(group.vec))$p.value)})
            apply(samp.dat[,genes.use,drop=F],2,function(x) t.test(x ~ as.factor(group.vec)))
        },
        error = function(cond) {
            print(i)
            print(unique(group.vec))
            
            message(cond)
            message(sprintf('Only %d entries for cluster %d out of %d cells',sum(cluster.vec==i),i,nrow(samp.dat)))

            break
        })
        ## cluster.markers[[i]] <- data.frame(p.val,row.names=colnames(counts.mat))
        
        names(gene.tests) <- genes.use

        pval.res[[sprintf('%d %d',i,j)]] <- data.frame(pvalue = sapply(gene.tests, function(x) x$p.value),statistic = sapply(gene.tests,function(x) x$statistic))
            ## c(gene.tests$p.value,gene.tests$statistic)

        print(sprintf('Completed cluster %d %d',i,j))

        saveRDS(pval.res,'results/recount/markers/anatomy_ttest_markers_pairwise.RDS')
    }
}


## saveRDS(cluster.markers,'results/allen/markers/kmeans_wilcox_markers.RDS')
## saveRDS(pval.mat,'results/allen/markers/louvain_wilcox_markers.RDS')
## 
## clusters <- clusters[rownames(allen.dat)]
##
## cluster.dat <- list()
##
## for(i in 1:length(unique(clusters))){
##
##     group.vec <- rep('Outside',nrow(allen.dat))
##     group.vec[clusters==i] <- 'Inside'
##
##     p.val <- tryCatch(
##     {
##         sapply(1:ncol(allen.dat),function(x) {return(wilcox.test(allen.dat[,x] ~ as.factor(group.vec))$p.value)})
##      
##     },
##     error = function(cond) {
##         message(cond)
##         message(sprintf('Only %d entries for cluster %d out of %d cells',sum(group.vec[clusters==i]),i,nrow(allen.dat)))
##
##         break
##     })
##  
##     cluster.dat[[i]] <- data.frame(p.val,row.names=colnames(allen.dat))
## }
##  
## saveRDS(cluster.dat,'results/allen/markers/kmeans_8_wilcox.RDS')
