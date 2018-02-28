#!/bin/bash

#PBS -N dl_allen_data
#PBS -j oe
#PBS -m abe
#PBS -M spudman.jobs@gmail.com
#PBS -o /gpfs/home/jbrugg/map/reuslts/logs/

WORK_DIR=/gpfs/group/su/lhgioia/map/

while read out_file get_url
do

    # if [ ! -f ${WORK_DIR}data/allen/gene_count_files/${out_file}.tsv ]
    if [[ $(wc -l <${WORK_DIR}data/allen/gene_count_files/${out_file}.tsv) -le 40000 ]]
    then
	echo http://api.brain-map.org${get_url}
	curl -o ${WORK_DIR}data/allen/gene_count_files/$out_file.tsv http://api.brain-map.org${get_url}
    fi

done < ${WORK_DIR}data/allen/tmp/allen_dl_info.tsv
