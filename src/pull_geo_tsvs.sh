#!/bin/bash

WORK_DIR='/gpfs/group/su/lhgioia/map/'

cd $PBS_O_WORKDIR

while read line
do

    id=${line%%_*}
    ftp_line=`grep $id data/geo_searches/ftp_list.csv`
    # echo $ftp_line

    curl ${ftp_line#*,}suppl/$line -o ${WORK_DIR}data/geo_searches/encode_tsvs/${line}

done < ${WORK_DIR}data/geo_searches/geo_tsvs_encode_search.txt

