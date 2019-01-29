#!/bin/bash

#PBS -N allen_centers
#PBS -q workq
#PBS -l mem=31gb
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -M spudman.jobs@gmail.com
#PBS -m abe

cd $PBS_O_WORKDIR

module load R/3.4.0

Rscript src/tmp/get_allen_centers.R

echo 'done'
