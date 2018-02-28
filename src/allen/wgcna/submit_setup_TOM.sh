#!/bin/bash

#PBS -N setup_TOM_rotation_genes_5k
#PBS -q bigmem
#PBS -l mem=251gb
#PBS -l walltime=400:00:00
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

module load R/3.4.3
module load gcc/6.3.0

cd $PBS_O_WORKDIR

Rscript src/allen/wgcna/setup_TOM.R
