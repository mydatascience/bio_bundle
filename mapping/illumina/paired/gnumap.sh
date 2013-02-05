#!/bin/bash
#USAGE ./gnumap.sh reference ref_base_name reads reads_base_name result threads hashsz additional

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
gnumap-plain $additional -g $ref -c $threads -m $hashsz -o $out/$reads1_base $reads1 $reads2 && \
echo "mapping started" && \
date
