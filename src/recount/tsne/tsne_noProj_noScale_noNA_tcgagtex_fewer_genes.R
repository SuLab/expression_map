library(Rtsne)

setwd('/gpfs/group/su/lhgioia/map')

## perplexity.vec <- c(10,20,30,50,75,100,500)
perplexity.vec <- c(30,40,60,80,110,140,175,210,250)

dim.vec <- c(50,100,150,200,250)

recount.data <- readRDS('results/recount/pca/recount_250_dim_noScaled_noProj_tcgagtex_fewer_genes.RDS')

## NA.rows <- read.table('data/recount/metadata/runs_not_in_metasra.csv',sep=',',header=F,stringsAsFactors=F)

## not.NAs <- rownames(recount.data)[!rownames(recount.data) %in% NA.rows$V1]

## recount.data <- recount.data[not.NAs,]

print(dim(recount.data))

results <- list()

for(perp in perplexity.vec){
    for(dim in dim.vec){
        tsne.out <- Rtsne(recount.data[,1:dim],perplexity=perp,check_duplicates=FALSE,pca=FALSE)
        rownames(tsne.out$Y) <- rownames(recount.data)
        results[[sprintf('%d.%d',perp,dim)]] <- tsne.out$Y

        saveRDS(results,'results/recount/tsne/recount_tsne_pca_list_noNAs_tcgagtex_fewer_genes.RDS')
        
    }
}
