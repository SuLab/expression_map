#!/bin/bash

cd /Users/Jake/Documents/Projects/Mercator/downloads/recount/

for file in `ls raw_counts/*`
do

    tail -n +2  $file | awk '{ for(i=1; i<=NF; i++) j+=$i; print j; j=0 }' | paste gene_counts.tsv - > foo
    mv foo gene_counts.tsv
    
done

