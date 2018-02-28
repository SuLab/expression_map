#!/bin/bash

#PBS -N remove_outlier_samples
#PBS -q bigmem
#PBS -l mem=251gb
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

module load R/3.4.0
module load gcc/6.3.0

cd $PBS_O_WORKDIR

Rscript src/allen/wgcna/remove_outlier_samples.R
