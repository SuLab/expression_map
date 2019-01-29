#!/bin/bash

#PBS -N kmeans_recount_sub
#PBS -q workq
#PBS -l mem=47gb
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R

Rscript src/recount/clustering/small_cluster_kmeans_recount.R
