#!/bin/bash

cd ~/map/

for file in `ls data/baldwin_data/adapter_clipped_reads/`
do

    rm -r results/allen/baldwin_rsem_first/${file%.fastq.gz}/

    mkdir results/allen/baldwin_rsem_first/${file%.fastq.gz}/

    # gunzip -c data/baldwin_data/adapter_clipped_reads/$file >> data/baldwin_data/unzipped_data/${file%.gz}

    qsub -N rsem_bowtie_${file%.fastq.gz} \
	 -v var="$file" \
	 src/allen/rsem_first.sh

	 # rsem-calculate-expression -p 8 \
	 # --bowtie-m 100 \
	 # --bowtie-e 500 \
	 # data/baldwin_data/unzipped_data/${file%.gz} \
	 # data/allen/rsem_gtf/ref/rsem_GRCm38 \
	 # results/allen/baldwin_rsem_first/${file%.fastq.gz}/${file%.fastq.gz}
    
done

	 # -q workq \
	 # -l nodes=1:ppn=8 \
	 # -l mem=47gb \
	 # -l walltime=50:00:00 \
	 # -o /gpfs/home/jbrugg/map/results/logs/ \
	 # -j oe \
	 # -m abe \
	 # -M spudman.jobs@gmail.com \
