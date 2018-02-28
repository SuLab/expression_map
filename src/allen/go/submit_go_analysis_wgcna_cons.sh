#!/bin/bash

#PBS -N go_analysis_wilvain
#PBS -q workq
#PBS -l mem=15gb
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.0
module load gcc/6.3.0

Rscript src/allen/go/go_analysis_wilvain_pairwise.R
