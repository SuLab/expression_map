#!/bin/sh

#PBS -N norm_samples_gtex
#PBS -q workq
#PBS -l nodes=1:ppn=1
#PBS -l mem=3gb
#PBS -l walltime=100:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

module load R/3.4.0

Rscript /gpfs/home/jbrugg/map/src/norm_samples.R
