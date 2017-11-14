#!/bin/bash
## create tcga phenotype table with only relevant annotations

cd /gpfs/group/su/lhgioia/map/data/recount

## create header
sed 1,3d pheno_col_matching.txt | cut -f 2 | paste -sd '\t' > tcga/pheno_relevant_header.tsv

cd /gpfs/group/su/lhgioia/map/data/recount/tcga

## create table body
awk -v FS='\t' '{$2=$1; $1=$19; gsub(".bw", ""); $9=$108; gsub("Solid Tissue ", ""); $4="NA"; $16="NA"; $19="NA"; $20="NA"; print $1, $2, $65, $4, $77, $74, $78, $90, $9, $92, $91, $67, $70, $40, $86, $16, $29, $193, $19, $20}' OFS='\t' phenotype.tsv | sed 1d > pheno_relevant_body.tsv

## cat header and body together
cat pheno_relevant_header.tsv pheno_relevant_body.tsv > pheno_relevant.tsv

## reset tumor info for normal samples
awk -v FS='\t' '{if($9=="Normal"){$6=$5; $7="NA"; $8="NA"; print $0;}else print $0}' OFS='\t' pheno_relevant.tsv > pheno_formatted.tsv


exit=0
