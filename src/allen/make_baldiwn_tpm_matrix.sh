#!/bin/bash

#PBS -N make_baldwin_tpm_matrix
#PBS -j oe
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

# cd $PBS_O_WORKDIR

cut -f 1 results/allen/baldwin_rsem_first/A1-B2-3/A1-B2-3.genes.results | tail -n +2 | tr '\r\n' ',' > data/allen/baldwin_tpm_matrix.tmp
echo '' >> data/allen/baldwin_tpm_matrix.tmp

for file in `ls results/allen/baldwin_rsem_first/`
do

    echo -n ${file}',' >> data/allen/baldwin_tpm_matrix.tmp
    cut -f 6 results/allen/baldwin_rsem_first/${file}/${file}.genes.results | tail -n +2 | tr '\r\n' ',' >> data/allen/baldwin_tpm_matrix.tmp
    echo '' >> data/allen/baldwin_tpm_matrix.tmp

done

sed -e 's/,$//' data/allen/baldwin_tpm_matrix.tmp > data/allen/baldwin_tpm_matrix.csv
