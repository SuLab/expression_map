#!/bin/bash

#PBS -N marker_gene_pval_table
#PBS -q workq
#PBS -l mem=7gbgb
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R

Rscript src/allen/sql/marker_gene_pval_table.R

echo 'done'
