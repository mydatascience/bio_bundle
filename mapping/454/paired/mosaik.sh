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

dir=$(dirname $(readlink -e $0))

echo "build started" 
date 
echo -e "0\tBuilding index\t"`date +%s` >> $out/mapping_time.log
MosaikBuild -fr $ref -oa $out/$ref_base.dat && \
MosaikBuild -q $reads1 -out $out/$reads1_base.dat -st 454 && \
MosaikBuild -q $reads2 -out $out/$reads2_base.dat -st 454 && \
echo "align started" && \
date && \
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
MosaikAligner -p $threads -hs $hashsz -in $out/$reads1_base.dat -in $out/$reads2_base.dat -out $out/$reads1_base -ia $out/$ref_base.dat -annpe $dir/../../../misc/pe.ann -annse $dir/../../../misc/se.ann && \
echo "align done" && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log && \
date
