#!/bin/bash
#USAGE ./gnumap.sh reference ref_base_name reads1 reads1_base_name reads2 reads2_base_name result threads hashsz additional

ref=$1
ref_base=$2
reads1=$3
reads1_base=$4
reads2=$5
reads2_base=$6
out=$7
threads=$8
hashsz=$9
additional=${10}

echo "mapping started" 
date
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
gnumap-plain $additional -g $ref -c $threads -m $hashsz -o $out/$reads1_base $reads1 $reads2 && \
echo "mapping started" && \
date &&
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log
