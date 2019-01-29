wilcox.markers <- readRDS('results/allen/markers/louvain_wilcox_markers.RDS')
zero.prop <- readRDS('results/allen/gene_zero_entries_cluster.RDS')

## gene.names <- read.table('data/allen/tmp/gene_list.txt')
## gene.names <- gsub('$','X',gene.names$V1)

## pval.df <- NULL
## stat.df <- NULL

all.markers <- NULL

marker.list <- list()

median.fc <- readRDS('results/allen/markers/louvain_median_fold_change_pairwise.RDS')
mean.fc <- readRDS('results/allen/markers/louvain_mean_fold_change_pairwise.RDS')

var.genes <- readRDS('data/allen/tmp/variable_genes_allen.RDS')

pairwise.cluster.pvals <- readRDS('results/allen/markers/louvain_wilcox_pairwise_compared_markers_all_genes.RDS')

for(i in 1:length(wilcox.markers)){

    pval.col <- stat.col <- rep(NA,length(zero.prop[[1]]))
    names(pval.col) <- names(stat.col) <- names(zero.prop[[1]])

    pvals <- p.adjust(wilcox.markers[[i]][1,],method='bonferroni')
    markers.pvals <- wilcox.markers[[i]][1,][pvals < 0.01]
    markers.stats <- wilcox.markers[[i]][2,][pvals < 0.01]

    zero.vec <- zero.prop[[i]][names(markers.pvals)]

    pairwise.pvals <- p.adjust(pairwise.cluster.pvals[[i]])[names(markers.pvals)]

    ## if(is.null(all.markers)){
    ##     all.markers <- data.frame(pval=markers.pvals,gene=names(markers.pvals),cluster=rep(i,length(markers.pvals)),zero.prop=zero.vec,stat=markers.stats)
    ## }else{
    ##     all.markers <- rbind(all.markers,data.frame(pval=markers.pvals,gene=names(markers.pvals),cluster=rep(i,length(markers.pvals)),zero.prop=zero.vec,stat=markers.stats))
    ## }

    if(is.null(all.markers)){
        all.markers <- data.frame(pval=markers.pvals,gene=names(markers.pvals),cluster=rep(i,length(markers.pvals)),zero.prop=zero.vec,stat=markers.stats,pairwise.pvals=pairwise.pvals)
    }else{
        all.markers <- rbind(all.markers,data.frame(pval=markers.pvals,gene=names(markers.pvals),cluster=rep(i,length(markers.pvals)),zero.prop=zero.vec,stat=markers.stats,pairwise.pvals=pairwise.pvals))
    }

    marker.list[[i]] <- wilcox.markers[[i]][,pvals < 0.01]

    ## pval.col[names(wilcox.markers[[i]])] <- wilcox.markers[[i]][1,]
    ## stat.col[names(wilcox.markers[[i]])] <- wilcox.markers[[i]][2,]

    ## if(is.null(marker.df)){
    ##     pval.df <- 
        
    ## }
    ## else{
    ## }
}

saveRDS(all.markers,'results/allen/markers/all_markers_all_genes_onevsall.RDS')
saveRDS(marker.list,'results/allen/markers/marker_list_all_genes_onevsall.RDS')

median.col <- rep(NA,nrow(all.markers))
median.abs.col <- rep(NA,nrow(all.markers))
mean.col <- rep(NA,nrow(all.markers))
mean.abs.col <- rep(NA,nrow(all.markers))
## median.col <- NULL

for(i in 1:nrow(all.markers)){

    row <- all.markers[i,]

    median.fc.vec <- c()
    mean.fc.vec <- c()

    for(j in 1:length(wilcox.markers)){

        if(row$cluster != j){

            median.fc.entry <- median.fc[[sprintf('%d %d',row$cluster,j)]]
            mean.fc.entry <- mean.fc[[sprintf('%d %d',row$cluster,j)]]
            
            names(median.fc.entry) <- var.genes
            names(mean.fc.entry) <- var.genes
            ## names(fc.entry) <- gsub('^','X',var.genes)

            median.fc.vec <- c(median.fc.vec,median.fc.entry[as.character(row$gene)])
            mean.fc.vec <- c(mean.fc.vec,mean.fc.entry[as.character(row$gene)])

        }

    }

    ## which.min(abs(fc.vec))

    ## if(is.null(median.col)){
    ##     median.col <- data.frame(fc.vec[which.min(abs(fc.vec))],i)
    ## }else{
    ##     median.col <- rbind(median.col,data.frame(fc.vec[which.min(abs(fc.vec))],i))
    ## }

    ## print(sprintf('row %d :',i))
    ## print(fc.vec[which.min(abs(fc.vec))])
    ## print(row)
    ## print('')

    median.col[i] <- median.fc.vec[which.min(median.fc.vec)]
    mean.col[i] <- mean.fc.vec[which.min(mean.fc.vec)]
    median.abs.col[i] <- median.fc.vec[which.min(abs(median.fc.vec))]
    mean.abs.col[i] <- mean.fc.vec[which.min(abs(mean.fc.vec))]

    ## median.col <- c(median.col,fc.vec[which.min(abs(fc.vec))])

    ## if((i %% 1000) == 0){
    ##     print(i)
    ## }



}

fc.vecs <- list()

fc.vecs['median.fc'] <- median.col
fc.vecs['mean.fc'] <- mean.col
fc.vecs['median.abs.col'] <- median.abs.col
fc.vecs['mean.abs.col'] <- mean.abs.col

saveRDS(fc.vecs,'results/allen/markers/all_markers_min_fc_change.RDS')

all.markers$median.fc.change <- median.col
all.markers$median.abs.fc.change <- median.abs.col
all.markers$mean.fc.change <- mean.col
all.markers$mean.abs.fc.change <- mean.abs.col

saveRDS(all.markers,'results/allen/markers/all_markers_all_genes_onevsall.RDS')




