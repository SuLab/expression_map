#!/bin/bash

#PBS -N welchanova_louvain_markers_all_genes
#PBS -q workq
#PBS -l mem=47gb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.3
module load gcc/6.3.0

Rscript src/allen/markers/kruskal_louvain_markers_all_genes.R
