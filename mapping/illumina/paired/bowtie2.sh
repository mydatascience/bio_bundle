#!/bin/bash
#USAGE ./bowtie2.sh reference ref_base_name reads1 reads1_base_name reads2 reads2_base_name result threads hashsz additional

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

echo "build started" 
date 
bowtie2-build $ref $out/$ref_base.ind && \
echo "mapping started" && \
date && \
bowtie2 -t -q -p $threads -a --fr --minins 0 --maxins 500 $out/$ref_base.ind -1 $reads1 -2 $reads2 > $out/$reads1_base.sam && \
echo "mapping done" && \
date
