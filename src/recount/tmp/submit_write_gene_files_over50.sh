#!/bin/bash

#PBS -N write_gene_files
#PBS -q workq
#PBS -l mem=47gb
#PBS -l walltime=100:00:00
#PBS -m abe
#PBS -M spudman.jobs@gmail.com
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/

cd $PBS_O_WORKDIR

module load R

Rscript src/recount/tmp/write_gene_files_over50.R
