#!/bin/bash

#PBS -N write_recount_tpm
#PBS -q workq
#PBS -l mem=47gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R

Rscript src/recount/tmp/write_tpm_matrix.R
