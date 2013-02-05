#!/bin/bash
#USAGE ./razers.sh reference ref_base_name reads reads_base_name result threads hashsz additional

ref=$1
ref_base=$2
reads=$3
reads_base=$4
out=$5
threads=$6
hashsz=$7
additional=$8

echo "mapping started" && \
date && \
razers3 $additional -i 94 -rr 97 -tc $threads --output-format sam -o $out/$reads_base.sam $ref $reads && \
echo "mapping done" && \
date
