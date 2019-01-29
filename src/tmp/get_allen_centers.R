setwd('/gpfs/group/su/lhgioia/map/')

allen.dat <- read.table('data/allen/tpm_matrix.csv',sep=',',header=T)

col.centers <- apply(allen.dat,2,mean)

saveRDS(col.centers,'data/allen/allen_center.RDS')
