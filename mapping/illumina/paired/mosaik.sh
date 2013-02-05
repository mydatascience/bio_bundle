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

echo "build started" 
date 
MosaikBuild -fr $ref -oa $out/$ref_base.dat && \
MosaikBuild -q $reads1 -out $out/$reads1_base.dat -st Illumina && \
MosaikBuild -q $reads2 -out $out/$reads2_base.dat -st Illumina && \
echo "align started" && \
date && \
MosaikAligner -p $threads -hs $hashsz -in $out/$reads1_base.dat -in $out/$reads2_base.dat -out $out/$reads1_base -ia $out/$ref_base.dat -annpe $MOSAIK_ANN/pe.ann -annse $MOSAIK_ANN/se.ann && \
echo "align done" && \
date
