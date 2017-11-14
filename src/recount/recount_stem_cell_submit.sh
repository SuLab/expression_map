#!/bin/bash

#PBS -N stem_cell_dl
#PBS -m abe
#PBS -M spudman.jobs@gmail.com
#PBS -o /gpfs/home/jbrugg/map/src/logs/

cd $PBS_O_WORKDIR

Rscript src/recount/recount_stem_cell.R
