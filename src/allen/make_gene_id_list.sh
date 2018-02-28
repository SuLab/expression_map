#!/bin/bash

# THI DOESNT WORK BC $line matches to the start/end positions of the genes, aka not returning the proper genes. need to change to matching 'gene_id "$line";

cd /gpfs/group/su/lhgioia/map

echo gene_id,transcript_id,gene_symbol > data/allen/gene_id_table.tmp

for line in `cat data/allen/gene_list.txt`
do
    grep -m 1 "gene_id \"${line}\";" data/allen/rsem_gtf/rsem_transfirst_GRCm38.p3.gtf | cut -f 9 >> data/allen/gene_id_table.tmp
done

sed -e 's/; /,/g' -e 's/gene_id //' -e 's/transcript_id //' -e 's/gene_symbol //' -e 's/"//g' data/allen/gene_id_table.tmp > data/allen/gene_id_table.csv
