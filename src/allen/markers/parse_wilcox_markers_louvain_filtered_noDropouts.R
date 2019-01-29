setwd('/gpfs/group/su/lhgioia/map/')

marker.dat <- readRDS('results/allen/markers/louvain_wilcox_markers_pairwise_noDropouts_filtered.RDS')

## gene.table <- read.table('data/allen/gene_id_table.csv',sep=',',header=T)
gene.ids <- read.table('data/allen/gene_list.txt')
gene.ids <- sapply(gene.ids$V1,function(x) sprintf('X%d',x))
## gene.ids <- sapply(gene.table$gene_id,function(x) sprintf('X%d',x))

## gene.ids <- readRDS('results/allen/variable_genes_log_4k.RDS')

cluster.pvals <- list()

for(i in 1:30){

    print(i)

    ## comp.mat <- matrix(0,nrow=nrow(gene.table),ncol=18)
    comp.mat <- matrix(1,nrow=length(gene.ids),ncol=29)

    cluster.comps <- names(marker.dat)[grepl(sprintf('^%d ',i),names(marker.dat)) | grepl(sprintf(' %d$',i),names(marker.dat))]

    if(length(cluster.comps) != 29){
        print(sprintf('Wrong number of comps for cluster %d',i))
        break
    }

    for(j in 1:length(cluster.comps)){
        comp.vec <- marker.dat[[cluster.comps[j]]][gene.ids,,drop=F]
        rownames(comp.vec) <- gene.ids
        comp.mat[,j] <- comp.vec$pvalue
    }
    
    rownames(comp.mat) <- gene.ids

    gene.pvals <- apply(comp.mat,1,max)

    cluster.pvals[[i]] <- gene.pvals

    saveRDS(cluster.pvals,'results/allen/markers/louvain_wilcox_pairwise_compared_markers_filtered_noDropouts.RDS')

}

## wilcox.markers <- readRDS('results/allen/markers/louvain_wilcox_pairwise_compared_markers.RDS')
gene.names <- read.table('data/allen/gene_list.txt')$V1
## gene.names <- readRDS('results/allen/variable_genes_log_4k.RDS')

marker.list <- list()

## gene.desigs <- as.list(rep(NA,nrow(gene.names)))
## gene.desigs <- rep(NA,nrow(gene.names))
gene.desigs <- rep(NA,length(gene.names))
## names(gene.desigs) <- gsub('^','X',gene.names$V1)
names(gene.desigs) <- gsub('^','X',gene.names)

for(i in 1:length(cluster.pvals)){
    print(i)
    pvals <- p.adjust(cluster.pvals[[i]],method='bonferroni')
    ## marker.list[[i]] <- pvals
    markers <- pvals[!is.na(pvals)]
    markers <- markers[markers < 0.01 ]
    markers <- sort(markers)

    marker.list[[i]] <- markers

    for(gene in names(markers)){

        if(is.na(gene.desigs[gene])){
            gene.desigs[gene] <- i
        }else{
            gene.desigs[gene] <- paste(gene.desigs[gene],i,sep=',')
        }
    }
}

saveRDS(gene.desigs,'results/allen/markers/pairwise_wilvain_gene_designations_filtered_noDropouts.RDS')
saveRDS(marker.list,'results/allen/markers/pairwise_wilvain_marker_list_filtered_noDropouts.RDS')

