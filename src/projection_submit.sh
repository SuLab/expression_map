#!/bin/bash

#PBS -N proj_big
#PBS -q bigmem
#PBS -l nodes=1:ppn=1
#PBS -l mem=252gb
#PBS -l walltime=100:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load python

python src/random_projection.py
