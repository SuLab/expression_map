library(ggplot2)

setwd('/gpfs/group/su/lhgioia/map')

## cluster.markers <- readRDS('results/allen/markers/kmeans_8_wilcox.RDS')

## filtered.markers <- list()

## for(i in 1:length(cluster.markers)){

##     markers <- cluster.markers[[i]][!is.na(cluster.markers[[i]]$p.val),]
##     names(markers) <- rownames(cluster.markers[[i]])[!is.na(cluster.markers[[i]]$p.val)]

##     markers <- markers[order(markers)]

##     threshold <- 0.01/length(markers)

##     ## passed.genes <- rownames(cluster.markers[[i]])[cluster.markers[[i]]$p.val < threshold]

##     passed.genes <- names(markers)[markers < threshold]

##     passed.genes <- gsub('X','',passed.genes)

##     filtered.markers[[i]] <- passed.genes

## }

## saveRDS(filtered.markers,'results/allen/markers/kmeans_10_marker_genes.RDS')

filtered.markers <- readRDS('results/allen/markers/kmeans_10_marker_genes.RDS')

tsne.dat <- readRDS('results/allen/tsne/allen_tsne.RDS')
tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)

for(i in 1:length(filtered.markers)){

    markers <- filtered.markers[[i]]

    for(marker in markers[1:10]){

        gene.vals <- tpm.dat[,sprintf('X%s',marker)]
        
        plot.dat <- data.frame(tsne1=tsne.dat$Y[,1],tsne2=tsne.dat$Y[,2],gene.val=gene.vals)

        g <- ggplot() +
            geom_point(data=plot.dat,aes(x=tsne1,y=tsne2,colour=gene.val)) +
            ggtitle(marker)

        ggsave(sprintf('results/plots/kmeans_marker_plots/kmean_marker_%s_%d.png',marker,i),g)

    }

}
