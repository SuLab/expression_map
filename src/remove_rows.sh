#!/bin/bash
## NOT WORTH IT!!! Only removes 130 lines from first 100 gtex samples
## remove rows with all zeros
file=""

awk '{for(i=1;i<=NF;i++) t+=$i; if(t!=0) print $0; t=0}' ${file}.tsv > ${file}_no0.tsv

exit=0
