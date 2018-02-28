#!/bin/bash

#PBS -N project_eigengenes_pheno_sep
#PBS -q workq
#PBS -l mem=47gb
#PBS -l walltime=400:00:00
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

module load R/3.4.0
module load gcc/6.3.0

cd $PBS_O_WORKDIR

Rscript src/allen/wgcna/project_eigengenes_separate_pheno_filtered.R
