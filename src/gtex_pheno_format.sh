#!/bin/bash
## create tcga phenotype table with only relevant annotations

cd /gpfs/group/su/lhgioia/map/data/recount

## create header
sed 1,3d pheno_col_matching.txt | cut -f 2 | paste -sd '\t' > gtex/pheno_relevant_header.tsv

cd /gpfs/group/su/lhgioia/map/data/recount/gtex

## create table body
awk -v FS='\t' '{$1=$4; $4=$2; $2="gtex"; $7="NA"; $8="NA"; $9="Normal"; $10="dead"; $11="NA"; $12="NA"; $13="NA"; $17="NA"; $18="NA"; print $1, $2, $20, $4, $25, $26, $7, $8, $9, $10, $11, $12, $13, $22, $22, $31, $17, $18, $21, $23}' OFS='\t' phenotype.tsv | sed 1d > pheno_relevant_body.tsv

## cat header and body together
cat pheno_relevant_header.tsv pheno_relevant_body.tsv > pheno_relevant.tsv

## add NA to last column if empty
awk -v FS='\t' '{if($20==""){$20="NA";} print $0}' OFS='\t' pheno_relevant.tsv > pheno_formatted.tsv

## create sample key for sex and age merging
cut -f3 pheno_formatted.tsv | cut -d "-" -f1,2 | paste pheno_formatted.tsv - > pheno_formatted_withKey.tsv

## change sex values in sex&age table
cut -f 1-3 GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt | awk -v FS='\t' '{if($2==1) {$2="male"}else $2="female"; print $0}' OFS='\t' > gtex_sex_age.tsv

## merge in R
module load R
Rscript /gpfs/group/su/lhgioia/map/src/gtex_pheno_merge.R


exit=0
