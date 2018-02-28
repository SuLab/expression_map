#!/bin/bash

#PBS -N kmeans_markers_allen
#PBS -q workq
#PBS -l mem=31gb
#PBS -l walltime=100:00:00
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.0

Rscript src/allen/markers/wilcox_kmeans_markers.R
