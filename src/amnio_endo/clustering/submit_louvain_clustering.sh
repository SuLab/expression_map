#!/bin/bash

#PBS -N louvain_clustering_recount_noAmnio_k30
#PBS -q ssd
#PBS -l mem=250gb
#PBS -l walltime=150:00:00
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.3
module load gcc/6.3.0

Rscript src/amnio_endo/clustering/louvain_clustering.R
