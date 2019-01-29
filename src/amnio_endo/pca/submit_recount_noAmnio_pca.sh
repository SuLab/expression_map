#!/bin/bash

#PBS -N recount_noAmnio_pca_scaled
#PBS -q workq
#PBS -l nodes=1:ppn=1
#PBS -l mem=47gb
#PBS -l walltime=10:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com


cd $PBS_O_WORKDIR

module load R
# module load R/3.4.3
module load gcc

Rscript src/amnio_endo/pca/recount_noAmnio_pca_scaled.R
