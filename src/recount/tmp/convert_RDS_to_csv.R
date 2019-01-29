recount.pca <- readRDS('results/recount/pca/recount_250_dim_noScaled_noProj_over50_noSingle.RDS')

write.table(recount.pca,'results/recount/pca/recount_250_dim_noScaled_noProj_over50_noSingle.csv',sep=',',row.names=F,col.names=F)
