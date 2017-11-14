#!/bin/bash

#PBS -N tsne_recount
#PBS -q workq
#PBS -l nodes=1:ppn=8
#PBS -l mem=31gb
#PBS -l walltime=100:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R

Rscript src/tsne_runs/recount_tsne.R
