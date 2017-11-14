#!/bin/bash

cd /gpfs/group/su/lhgioia/map/data/recount/pheno

## remove old results file
rm tcga_organ_in_metasra.tsv

## get unique organ entries from tcga data
grep "TCGA" all_recount_metasra_summarized.tsv | sed 1d | cut -f6 | sort | uniq > tcga_organ_uniq.txt

## run through unique tcga organ entries and grep from metadata file
filename='tcga_organ_uniq.txt'
while read p; do
    # echo $p 
    grep "$p" all_recount_metasra_summarized.tsv >> tcga_organ_in_metasra.tsv
done < $filename

grep -v "SRP012682" tcga_organ_in_metasra.tsv > tcga_organ_in_metasra_noGtex.tsv

grep "SRP012682" all_recount_metasra_summarized.tsv > gtex_in_metasra.tsv

grep "stem cell" all_recount_metasra_summarized.tsv > stemcell_in_metasra.tsv

cat all_recount_metasra_summarized.tsv.header gtex_in_metasra.tsv tcga_organ_in_metasra_noGtex.tsv stemcell_in_metasra.tsv > recount_metasra_gtex_tcga_stem.tsv
