#!/bin/bash

#PBS -N allen_euclid_dist
#PBS -q workq
#PBS -l mem=31gb
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.0

Rscript src/allen/euclid_distance_allen_baldwin.R

echo 'done'
