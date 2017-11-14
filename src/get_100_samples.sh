#!/bin/bash
## removes row names and takes first 100 samples
cd /gpfs/group/su/lhgioia/map/data/recount/scaled

for dataset in gtex tcga; do
	head -n1 ${dataset}_tpm.tsv > ${dataset}_header.tsv
	sed 1d ${dataset}_tpm.tsv | cut -d " " -f 2- > ${dataset}_body.tsv
	cat ${dataset}_header.tsv ${dataset}_body.tsv > ${dataset}_tpm_norn.tsv
	cut -d " " -f 1-100 ${dataset}_tpm_norn.tsv > ${dataset}_tpm_100.tsv
done


exit=0
