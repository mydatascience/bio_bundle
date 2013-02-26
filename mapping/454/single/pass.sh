#!/bin/bash
#USAGE ./pass.sh reference ref_base_name reads reads_base_name result threads hashsz additional

ref=$1
ref_base=$2
reads=$3
reads_base=$4
out=$5
threads=$6
hashsz=$7
additional=$8

echo "mapping started"
date
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
pass -p 1111110111111 $additional -cpu $threads -fid 90 -sam -phred33 -fastq $reads -d $ref > $out/$reads_base.sam && \
echo "mapping done" && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log && \
date
