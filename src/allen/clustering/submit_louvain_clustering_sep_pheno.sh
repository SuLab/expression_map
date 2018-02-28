#!/bin/bash

#PBS -N louvain_clustering_allen_filtered_k3
#PBS -q workq
#PBS -l mem=31gb
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.0
module load gcc/6.3.0

Rscript src/allen/clustering/louvain_clustering_sep_pheno.R
