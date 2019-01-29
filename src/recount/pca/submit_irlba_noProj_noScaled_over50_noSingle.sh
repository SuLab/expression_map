#!/bin/bash

#PBS -N recount_noProj_noSingle_pca
#PBS -q workq
#PBS -l nodes=1:ppn=1
#PBS -l mem=251gb
#PBS -l walltime=100:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R
module load gcc

Rscript src/recount/pca/irlba_noProj_noScaled_over50_noSingle.R



