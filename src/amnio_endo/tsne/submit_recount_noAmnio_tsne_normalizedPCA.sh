#!/bin/bash

#PBS -N noAnnio_tsne_scaled
#PBS -q workq
#PBS -l nodes=1:ppn=8
#PBS -l mem=31gb
#PBS -l walltime=100:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.3

# Rscript src/amnio_endo/tsne/recount_noNA_pca.R
Rscript src/amnio_endo/tsne/recount_noAmnio_tsne_normalizedPCA.R

