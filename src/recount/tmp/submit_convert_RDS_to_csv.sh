#!/bin/bash

#PBS -N convert_RDS_2_csv
#PBS -q workq
#PBS -l mem=47gb
#PBS -l walltime=100:00:00
#PBS -m abe
#PBS -M spudman.jobs@gmail.com
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/

cd $PBS_O_WORKDIR

module load R

Rscript src/recount/tmp/convert_RDS_to_csv.R
