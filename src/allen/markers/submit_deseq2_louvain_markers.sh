#!/bin/bash

#PBS -N deseq2_louvain_markers_allen_cluster1
#PBS -q workq
#PBS -l mem=125gb
#PBS -l nodes=1:ppn=16
#PBS -l walltime=400:00:00
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.3
module load gcc/6.3.0

Rscript src/allen/markers/deseq2_louvain_markers.R
