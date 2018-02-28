#!/bin/bash

#PBS -N make_allen_expected_matrix
#PBS -j oe
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

for gene_file in `ls data/allen/gene_count_files/`
do

    echo -n $gene_file',' >> data/allen/tpm_matrix.tsv
    cut -f 6 data/allen/gene_count_files/${gene_file} | tail -n +2 | tr '\r\n' ',' >> data/allen/tpm_matrix.tsv
    echo '' >> data/allen/tpm_matrix.tsv

done

