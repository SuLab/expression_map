#!/bin/bash

#PBS -N tsne_allen_filtered_log_4k
#PBS -q workq
#PBS -l nodes=1:ppn=8
#PBS -l mem=31gb
#PBS -l walltime=100:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R/3.4.0

Rscript src/allen/tsne/tsne_allen_filtered_log_4k.R
