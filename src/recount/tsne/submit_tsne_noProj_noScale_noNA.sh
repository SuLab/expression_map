#!/bin/bash

#PBS -N tsne_parameters_noNA
#PBS -q workq
#PBS -l nodes=1:ppn=8
#PBS -l mem=31gb
#PBS -l walltime=200:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load R

Rscript src/recount/tsne/tsne_noProj_noScale_noNA.R
