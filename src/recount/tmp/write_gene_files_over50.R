setwd('/gpfs/group/su/lhgioia/map')

tcga.tpm <- as.matrix(readRDS('data/recount/tcga/tpm_tcga.RDS'))
gtex.tpm <- as.matrix(readRDS('data/recount/gtex/tpm_gtex.RDS'))

## tpm.mat <- rbind(tcga.tpm,gtex.tpm)

tpm.mat <- matrix(0,nrow=70603,ncol=58037)

cnt <- 1

tpm.mat[cnt:nrow(tcga.tpm),] <- tcga.tpm
cnt <- cnt + nrow(tcga.tpm)

tpm.rownames <- rownames(tcga.tpm)

tpm.mat[cnt:(cnt + nrow(gtex.tpm)-1),] <- gtex.tpm
cnt <- cnt + nrow(gtex.tpm)

tpm.rownames <- c(tpm.rownames,rownames(gtex.tpm))

over50.projs <- readRDS('data/recount/recount_entries_over50_noSingle.RDS')

## for(file.id in list.files('data/recount/project_cnts')){
for(file.id in over50.projs){

    file.tpm <- as.matrix(readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_tpm.RDS',file.id)))

    tpm.mat[cnt:(cnt + nrow(file.tpm) - 1),] <- file.tpm
    cnt <- cnt + nrow(file.tpm)

    tpm.rownames <- c(tpm.rownames,rownames(file.tpm))
    ## tpm.mat <- rbind(tpm.mat,file.tpm)

}

tpm.mat <- tpm.mat[1:length(tpm.rownames),]

print(dim(tpm.mat))

rownames(tpm.mat) <- tpm.rownames

## colnames(tpm.mat) <- gsub('[.].*$','',colnames(tcga.tpm))
colnames(tpm.mat) <- colnames(tcga.tpm)

done.flag <- TRUE

i <- 1

while(done.flag){

    first.ind <- (i-1)*10000+1
    second.ind <- 10000*i

    if(second.ind > ncol(tpm.mat)){
        done.flag <- FALSE
        second.ind <- ncol(tpm.mat)
    }
        
    ## write.table(t(round(tpm.mat[first.ind:second.ind,],digits=6)),sprintf('data/recount/tmp/tpm_mat_over50_%d.csv',i),sep=',')
    write.table(t(round(tpm.mat[,first.ind:second.ind],digits=6)),sprintf('data/recount/tmp/tpm_mat_over50_%d.csv',i),sep=',')

    i <- i+1

}

saveRDS(rownames(tpm.mat),'data/recount/tmp/tpm_mat_rownames.RDS')
