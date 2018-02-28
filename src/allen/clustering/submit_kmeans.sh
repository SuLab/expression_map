#!/bin/bash

#PBS -N kmeans_allen_1to20_filtered_4k
#PBS -q bigmem
#PBS -l mem=251gb
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.0

Rscript src/allen/clustering/kmeans_allen.R
