library(WGCNA)

enableWGCNAThreads()

setwd('/gpfs/group/su/lhgioia/map')

options(stringsAsFactors=F)

tpm.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)
var.genes <- readRDS('results/allen/variable_genes_log_4k.RDS')
allen.meta <- read.table('data/allen/allen_meta.tsv',sep='\t',header=T)

tpm.dat <- tpm.dat[as.character(allen.meta$rnaseq_profile_id),var.genes]

mag.cnst <- as.integer(log(min(as.vector(tpm.dat)[as.vector(tpm.dat)!=0])))
dec.cnst <- exp(mag.cnst)
    
tpm.dat <- log(tpm.dat + dec.cnst) - mag.cnst

beta <- 2

expr.dat <- vector(mode='list',length=3)

cnt <- 1

for(i in unique(allen.meta$sources)){

    dat.slice <- tpm.dat[allen.meta$sources==i,]

    expr.dat[[cnt]] <- list(data=dat.slice)
    ## colnames(expr.dat[[cnt]]$data) <- colnames(dat.slice)
    ## rownames(expr.dat[[cnt]]$data) <- rownames(dat.slice)

    cnt <- cnt + 1
}

names(expr.dat) <- unique(allen.meta$sources)

exprSize <- checkSets(expr.dat)

if(!exprSize$structureOK){

    print('structure bad')
    quit()

}

adj.mat <- array(0,dim=c(3,exprSize$nGenes,exprSize$nGenes))

for(set in 1:3){
    
    adj.mat[set,,] <- abs(cor(expr.dat[[set]]$data,use='p'))^beta

}

TOM.all <- array(0,dim=c(3,exprSize$nGenes,exprSize$nGenes))

for(set in 1:3){

    TOM.all[set,,] <- TOMsimilarity(adj.mat[set,,])

}

scale.p <- 0.95
n.samples <- as.integer(1/(1-scale.p) * 1000)

scale.sample.ind <- sample(exprSize$nGenes*(exprSize$nGenes-1)/2,size=n.samples)

TOM.scaling.samples <- list()

scale.quant <- rep(1,3)
scale.powers <- rep(1,3)

for(i in 1:3){

    TOM.scaling.samples[[i]] <- as.dist(TOM.all[i,,])[scale.sample.ind]
    scale.quant[i] <- quantile(TOM.scaling.samples[[i]],probs=scale.p,type=8)

    if(i > 1){
        
        scale.powers[i] <- log(scale.quant[1])/log(scale.quant[i])
        TOM.all[i,,] <- TOM.all[i,,]^scale.powers[i]

    }

}

scaled.TOM.samples <- list()

for(i in 1:3){

    scaled.TOM.samples[[i]] <- TOM.scaling.samples[[i]]^scale.powers[i]

}

saveRDS(scaled.TOM.samples,'data/allen/tmp/scaled_tom_samples.RDS')
saveRDS(TOM.scaling.samples,'data/allen/tmp/unscaled_tom_samples.RDS')
saveRDS(TOM.all,'data/allen/wgcna/scaled_TOM_pheno.RDS')
