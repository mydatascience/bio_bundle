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
MosaikBuild -fr $ref -oa $out/$ref_base.dat && \
MosaikBuild -q $reads -out $out/$reads_base.dat -st 454 && \
echo "align started" && \
date && \
MosaikAligner -p $threads -hs $hashsz -in $out/$reads_base.dat -out $out/$reads_base -ia $out/$ref_base.dat -annpe $MOSAIK_ANN/pe.ann -annse $MOSAIK_ANN/se.ann && \
echo "align done" && \
date
