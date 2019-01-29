setwd('/gpfs/group/su/lhgioia/map/')

marker.dat <- readRDS('results/recount/markers/kmeans40_ttest_markers_pairwise.RDS')

## gene.table <- read.table('data/allen/gene_id_table.csv',sep=',',header=T)
## gene.ids <- read.table('data/allen/gene_list.txt')
## gene.ids <- sapply(gene.ids$V1,function(x) sprintf('X%d',x))
## gene.ids <- sapply(gene.table$gene_id,function(x) sprintf('X%d',x))

## gene.ids <- readRDS('results/allen/variable_genes_log_4k.RDS')

gene.file <- read.table('data/recount/gene_names_all.txt')
gene.ids <- gene.file$V1

cluster.pvals <- list()
cluster.fcs <- list()

for(i in 1:40){

    print(i)

    ## comp.mat <- matrix(0,nrow=nrow(gene.table),ncol=18)
    comp.mat <- matrix(1,nrow=length(gene.ids),ncol=39)
    fc.mat <- matrix(0,nrow=length(gene.ids),ncol=39)

    cluster.comps <- names(marker.dat)[grepl(sprintf('^%d ',i),names(marker.dat)) | grepl(sprintf(' %d$',i),names(marker.dat))]

    if(length(cluster.comps) != 39){
        print(sprintf('Wrong number of comps for cluster %d',i))
        break
    }

    for(j in 1:length(cluster.comps)){
        comp.vec <- marker.dat[[cluster.comps[j]]][gene.ids,,drop=F]
        rownames(comp.vec) <- gene.ids
        comp.mat[,j] <- comp.vec$pvalue
        if(i <= j){
            fc.mat[,j] <- -1 * comp.vec$fold.change
        }
        else{
            fc.mat[,j] <- comp.vec$fold.change
        }
        ## fc.mat[,j] <- comp.vec$fold.change
    }
    
    rownames(comp.mat) <- gene.ids
    rownames(fc.mat) <- gene.ids
    
    gene.pvals <- apply(comp.mat,1,max)
    gene.fcs <- apply(fc.mat,1,max)
    
    cluster.pvals[[i]] <- gene.pvals
    cluster.fcs[[i]] <- gene.fcs
    
    saveRDS(cluster.pvals,'results/recount/markers/kmeans_pairwise_compared_pvals.RDS')
    saveRDS(cluster.fcs,'results/recount/marekrs/kmeans_pairwise_compared_fcs.RDS')
}



## wilcox.markers <- readRDS('results/allen/markers/louvain_wilcox_pairwise_compared_markers.RDS')
## gene.names <- read.table('data/allen/gene_list.txt')$V1
## gene.names <- readRDS('results/allen/variable_genes_log_4k.RDS')

marker.list <- list()

## gene.desigs <- as.list(rep(NA,nrow(gene.names)))
## gene.desigs <- rep(NA,nrow(gene.names))
## gene.desigs <- rep(NA,length(gene.names))
## ## names(gene.desigs) <- gsub('^','X',gene.names$V1)
## names(gene.desigs) <- gsub('^','X',gene.names)

## gene.desigs <- gene.ids

for(i in 1:length(cluster.pvals)){
    print(i)
    pvals <- p.adjust(cluster.pvals[[i]],method='bonferroni')
    fcs <- cluster.fcs[[i]]
    ## marker.list[[i]] <- pvals
    markers <- pvals[!is.na(pvals)]
    fcs <- fcs[!is.na(pvals)]
    
    fcs <- fcs[markers < 0.01]
    markers <- markers[markers < 0.01 ]

    ## markers <- sort(markers)
    markers <- markers[order(fcs,decreasing=T)]

    marker.list[[i]] <- markers

    ## for(gene in names(markers)){

    ##     if(is.na(gene.desigs[gene])){
    ##         gene.desigs[gene] <- i
    ##     }else{
    ##         gene.desigs[gene] <- paste(gene.desigs[gene],i,sep=',')
    ##     }
    ## }
}

## saveRDS(gene.desigs,'results/recount/markers/pairwise_wilvain_gene_designations_filtered_noDropouts.RDS')
saveRDS(marker.list,'results/allen/markers/pairwise_wilvain_marker_list_filtered_noDropouts.RDS')

