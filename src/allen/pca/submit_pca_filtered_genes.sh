#!/bin/bash

#PBS -N pca_filtered_genes
#PBS -q workq
#PBS -l mem=31gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

module load R/3.4.0
module load gcc/6.3.0

cd $PBS_O_WORKDIR

Rscript src/allen/pca_filtered_genes.R
