#!/bin/bash

#PBS -N ttest_kmeans40_recount
#PBS -q workq
#PBS -l mem=93gb
#PBS -l walltime=150:00:00
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

module load R

cd $PBS_O_WORKDIR

Rscript src/recount/markers/kmeans_pairwise_markers.R

echo  'done'
