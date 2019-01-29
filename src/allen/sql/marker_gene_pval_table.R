setwd('/gpfs/group/su/lhgioia/map')

## marker.list <- readRDS('results/allen/markers/pairwise_wilvain_marker_list.RDS')
marker.list <- readRDS('results/allen/markers/all_markers_all_genes_onevsall.RDS')
gene.dat <- read.table('data/allen/gene_id_table.csv',sep=',',header=T)

gene.df <- data.frame()

for(i in 1:nrow(marker.list)){

    row <- marker.list[i,]

    gene <- gsub('X','',row$gene)
    gene.info <- gene.dat[gene.dat$gene_id==gene,]

    gene.entry <- data.frame(gene_id=gene.info$gene_id,
                             gene_symbol=gene.info$gene_symbol,
                             ## p_val=formatC(row$pairwise.pvals,format='e',digits=2),
                             p_val=row$pairwise.pvals,
                             cluster=row$cluster,
                             median_fc_change=row$median.fc.change,
                             median_abs_fc_change=row$median.abs.fc.change,
                             mean_fc_change=row$mean.fc.change,
                             mean_abs_fc_change=row$mean.fc.change,
                             ## one_vs_all_p_val=formatC(row$pval,format='e',digits=2),
                             one_vs_all_p_val=row$pval,
                             zero_prop=formatC(row$zero.prop,digits=3),
                             statistic=row$stat)


    if(nrow(gene.df) < 1){
        gene.df <- gene.entry
    }else{
        gene.df <- rbind(gene.df,gene.entry)
    }

}
    
write.table(gene.df,'results/allen/markers/gene_id_marker_table.csv',sep=',',row.names=F)


## order.fc <- function(i){

##     markers <- subset(all.markers,cluster==i)

##     ## head(markers[order(abs(markers$fc.change),decreasing=T),])

##     return(markers[order(markers$median.fc.change,decreasing=T),])
## }

## for(i in 1:length(marker.list)){
##     print(i)
##     j <- 0
##     for(gene in names(marker.list[[i]])){
##         j <- j + 1
##         print(j)
##         p.val <- marker.list[[i]][gene]        
##         gene <- gsub('X','',gene)
##         gene.info <- gene.dat[gene.dat$gene_id==gene,]
##         gene.entry <- data.frame(gene_id=gene.info$gene_id,gene_symbol=gene.info$gene_symbol,p_val=formatC(p.val,format='e',digits=2),cluster=i,rank=j)
##         if(nrow(gene.df) < 1){
##             gene.df <- gene.entry
##         }
##         else{
##             gene.df <- rbind(gene.df,gene.entry)
##         }
##     }
## }
