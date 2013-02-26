#!/bin/bash
#USAGE ./segemehl.sh reference ref_base_name reads reads_base_name result threads hashsz additional

ref=$1
ref_base=$2
reads=$3
reads_base=$4
out=$5
threads=$6
hashsz=$7
additional=$8

echo "build index started"
date 
echo -e "0\tBuilding index\t"`date +%s` >> $out/mapping_time.log
segemehl.x -t $threads -x $out/$ref_base.idx -d $ref && \
echo "mapping started" && \
date && \
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
segemehl.x $additional -t $threads -i $out/$ref_base.idx -d $ref -q $reads > $out/$reads_base.sam && \
echo "mapping done" && \
date && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log
