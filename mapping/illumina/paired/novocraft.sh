#!/bin/bash
#USAGE ./novocraft.sh reference ref_base_name reads1 reads1_base_name reads2 reads2_base_name result threads hashsz additional

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

echo "building index started" 
date 
novoindex -k $hashsz -t $threads $out/$ref_base $ref && \
echo "alignment started" && \
date && \
novoalign -d $out/$ref_base -f $reads1 $reads2 > $out/$reads1_base.sam && \
echo "alignment done" && \
date
