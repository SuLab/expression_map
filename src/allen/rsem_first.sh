#!/bin/bash

#PBS -q workq
#PBS -l nodes=1:ppn=8
#PBS -l walltime=100:00:00
#PBS -l mem=47gb
#PBS -o /gpfs/home/jbrugg/map/results/logs/
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com

cd $PBS_O_WORKDIR

module load bowtie/1.1.2

gunzip -c data/baldwin_data/adapter_clipped_reads/$var > data/baldwin_data/unzipped_data/${var%.gz}

# rsem-calculate-expression -p 8 \
# 			  --bowtie-m 100 \
# 			  --bowtie-e 500 \
# 			  data/baldwin_data/unzipped_data/Tac2-Hi-6.2.fastq \
# 			  data/allen/rsem_gtf/ref/rsem_GRCm38 \
# 			  results/allen/baldwin_rsem_first/Tac2-Hi-6.2/Tac2-Hi-6.2

rsem-calculate-expression -p 8 \
			  --bowtie-m 100 \
			  --bowtie-e 500 \
			  data/baldwin_data/unzipped_data/${var%.gz} \
			  data/allen/rsem_gtf/ref/rsem_GRCm38 \
			  results/allen/baldwin_rsem_first/${var%.fastq.gz}/${var%.fastq.gz}

			      # # --bowtie-path /opt/applications/bowtie/1.1.2/gnu/bin/bowtie \
			      # --gzipped-read-file \

