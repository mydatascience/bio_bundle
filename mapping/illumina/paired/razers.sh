#!/bin/bash
#USAGE ./mosaik.sh reference ref_base_name reads1 reads1_base_name reads2 reads2_base_name result threads hashsz additional

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

echo "mapping started" && \
date && \
razers3 $additional -i 94 -rr 97 -tc $threads --output-format sam -o $out/$reads1_base.sam $ref $reads1 $reads2 && \
echo "mapping done" && \
date
