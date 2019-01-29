#!/bin/bash

#PBS -N vptree_nn_recount
#PBS -q workq
#PBS -l mem=31gb
#PBS -l walltime=250:00:00
#PBS -j oe
#PBS -o /gpfs/group/su/lhgioia/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load python

python src/recount/clustering/nearest_neighbors_recount.py
