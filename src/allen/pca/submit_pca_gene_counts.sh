#!/bin/bash

#PBS -N nv500_pca_all_allen
#PBS -q workq
#PBS -l nodes=1:ppn=1
#PBS -l mem=47gb
#PBS -l walltime=150:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com


cd $PBS_O_WORKDIR

module load R/3.4.3

Rscript src/allen/pca/pca_gene_counts.R
# Rscript src/allen/pca_allen_baldwin.R
