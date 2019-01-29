#!/bin/bash

#PBS -N recount_noNAnoAmnio_pca.R
#PBS -q workq
#PBS -l nodes=1:ppn=1
#PBS -l mem=252gb
#PBS -l walltime=600:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com


cd $PBS_O_WORKDIR

module load R/3.4.3

Rscript src/amnio_endo/recount_noNA_pca.R
