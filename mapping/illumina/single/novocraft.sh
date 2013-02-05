#!/bin/bash
#USAGE ./novocraft.sh reference ref_base_name reads reads_base_name result threads hashsz additional

ref=$1
ref_base=$2
reads=$3
reads_base=$4
out=$5
threads=$6
hashsz=$7
additional=$8

echo "building index started" 
date 
novoindex -k $hashsz -t $threads $out/$ref_base $ref && \
echo "alignment started" && \
date && \
novoalign -d $out/$ref_base -f $reads > $out/$reads_base.sam && \
echo "alignment done" && \
date
