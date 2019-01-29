#!/bin/bash

#PBS -N kmeans_clustering_recount_noAmnio
#PBS -q workq
#PBS -l mem=250gb
#PBS -l walltime=150:00:00
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.3
module load gcc/6.3.0

Rscript src/amnio_endo/clustering/kmeans_clustering.R
