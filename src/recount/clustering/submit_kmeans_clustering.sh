#!/bin/bash

#PBS -N kmeans_recount
#PBS -q workq
#PBS -l mem=93gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R

Rscript src/recount/clustering/kmeans_recount.R
