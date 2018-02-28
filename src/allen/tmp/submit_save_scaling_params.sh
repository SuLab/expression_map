#!/bin/bash

#PBS -N save_scaling_params
#PBS -q bigmem
#PBS -l mem=251gb
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

module load R/3.4.0

cd $PBS_O_WORKDIR

Rscript src/allen/tmp/save_scaling_params.R
