#!/bin/bash
#USAGE ./bowtie.sh reference ref_base_name reads reads_base_name result threads hashsz additional

ref=$1
ref_base=$2
reads=$3
reads_base=$4
out=$5
threads=$6
hashsz=$7
additional=$8

date 
echo -e "0\tBuilding index\t"`date +%s` >> $out/mapping_time.log
tmap index -f $ref -w $hashsz && \
echo "mapping started" && \
date && \
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
tmap map3 $additional -R SM:Unknown -n $threads -f $ref -r $reads > $out/$reads_base.sam && \
echo "mapping done" && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log && \
date
