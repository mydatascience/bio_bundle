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

echo "mapping started" 
date 
tmap index -f $ref -w $hashsz && \
tmap map3 $additional -R SM:Unknown -n $threads -f $ref -r $reads > $out/$reads_base.sam && \
echo "mapping done" && \
date
