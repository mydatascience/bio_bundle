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

echo "build started" 
date 
bowtie-build $ref $out/$ref_base.ind && \
echo "mapping started" && \
date && \
bowtie -t -q -p $threads $additional --best --sam $out/$ref_base.ind $reads > $out/$reads_base.sam && \
echo "mapping done" && \
date
