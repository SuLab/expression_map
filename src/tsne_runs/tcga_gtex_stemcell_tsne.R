library(Rtsne)

setwd('/gpfs/group/su/lhgioia/map/')

stem.meta <- read.table('data/recount/pheno/recount_metasra_gtex_tcga_stem_uniq.tsv',sep='\t',header=T)

## rot.matrix <- as.matrix(read.table('data/recount/random_proj_rotation.csv',sep=','))

## print('rotation matrix read')

## tcga.tpm <- t(as.matrix(readRDS('data/recount/tcga/tpm_tcga.RDS')))
## tcga.rot <- t(rot.matrix %*% tcga.tpm)
## rownames(tcga.rot) <- colnames(tcga.tpm)
## saveRDS(tcga.rot,'data/recount/tcga/tcga_rot_tpm.RDS')

tcga.rot <- readRDS('data/recount/tcga/proj_tcga.RDS')

print('tcga read')

## gtex.tpm <- t(as.matrix(readRDS('data/recount/gtex/tpm_gtex.RDS')))
## gtex.rot <- t(rot.matrix %*% gtex.tpm)
## rownames(gtex.rot) <- colnames(gtex.tpm)
## saveRDS(gtex.rot,'data/recount/gtex/gtex_rot_tpm.RDS')

gtex.rot <- readRDS('data/recount/gtex/proj_gtex.RDS')

print('gtex read')

tsne.init <- rbind(gtex.rot,tcga.rot)

rm(gtex.rot)
rm(tcga.rot)
## rm(gtex.tpm)
## rm(tcga.tpm)
gc()

projects <- unique(stem.meta$proj_id)

projects <- projects[projects %in% list.files('data/recount/project_cnts/')]

for(proj.id in projects){
    
    print(proj.id)

    if(proj.id == 'TCGA' | proj.id == 'SRP012682'){
        next
    }

    ## proj.rot <- tryCatch({

    proj.rot <- readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_proj.RDS',proj.id))
    ##     return(proj.rot)
        
    ## }, error = function(err) {

    ##     print(sprintf('rotating %s',proj.id))

    ##     proj.tpm <- t(as.matrix(readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_tpm.RDS',proj.id))))
    ##     proj.rot <- t(rot.matrix %*% proj.tpm)
    ##     rownames(proj.rot) <- colnames(proj.tpm)
    ##     return(proj.rot)

    ## })
    
    tsne.init <- rbind(tsne.init,proj.rot)

}

## rm(rot.matrix)
## rm(proj.tpm)
## rm(proj.rot)
gc()

print('beginning tsne')

tsne.pca <- prcomp(tsne.init,retx=TRUE)
saveRDS(tsne.pca,'gpfs/home/jbrugg/map/results/recount/pca/tcga_gtex_stemcell_pca.RDS')

tsne.out <- Rtsne(tsne.pca$x[,1:50],perplexity=30, check_duplicates=FALSE,pca=FALSE)

saveRDS(tsne.out,file='/gpfs/home/jbrugg/map/results/recount/tsne/tsne_tcga_gtex_stemcell.RDS')
