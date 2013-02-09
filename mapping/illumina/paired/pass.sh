#!/bin/bash
#USAGE ./pass.sh reference ref_base_name reads1 reads1_base_name reads2 reads2_base_name result threads hashsz additional

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
pass -p 1111110111111 $additional -cpu $threads -fid 90 -block 50000 -sam -phred33 -fastq $reads1 -fastq $reads2 -d $ref > $out/$reads1_base.sam && \
echo "mapping done" && \
date
