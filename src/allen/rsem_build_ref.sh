#!/bin/bash

#PBS -N rsem_build_reference_chrnums
#PBS -q workq
#PBS -l nodes=1:ppn=1
#PBS -l mem=5gb
#PBS -l walltime=100:00:00
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load bowtie/1.1.2

rsem-prepare-reference --gtf data/allen/rsem_gtf/rsem_transfirst_chrnums_GRCm38.p3.gtf \
		       --bowtie \
		       doc/assemblies/mm10/GRCm38.p3.genome.fa \
		       data/allen/rsem_gtf/ref/rsem_GRCm38

echo 'done'

		       # --transcript-to-gene-map doc/assemblies/gene_transcript_trans.txt \

		       # --bowtie-path /opt/applications/bowtie/1.1.2/gnu/bin/ \
