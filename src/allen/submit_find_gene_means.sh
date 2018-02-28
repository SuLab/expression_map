#!/bin/bash

#PBS -N find_gene_means
#PBS -q workq
#PBS -l mem=31gb
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

module load R/3.4.0
module load gcc/6.3.0

cd $PBS_O_WORKDIR

Rscript src/allen/find_gene_means.R
