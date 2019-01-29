#!/bin/bash

#PBS -N euc_dist_recount
#PBS -q workq
#PBS -l nodes=1:ppn=16
#PBS -l mem=100gb
#PBS -l walltime=150:00:00
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.3
module load gcc/6.3.0

Rscript src/amnio_endo/map_project_samples_euclidian.R
