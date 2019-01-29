#!/bin/bash

#PBS -N read_trimming
#PBS -q workq
#PBS -l nodes=1:ppn=1
#PBS -l mem=6gb
#PBS -l walltime=400:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

for read_file in `ls /gpfs/home/jbrugg/baldwin_data/`
do

    ./lib/ea-utils/clipper/fastq-mcf doc/baldwin/adapters.fa /gpfs/home/jbrugg/baldwin_data/${read_file} -o data/baldwin_data/adapter_clipped_reads/${read_file}
done

