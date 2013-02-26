#!/bin/bash
#USAGE ./mosaik.sh reference ref_base_name reads reads_base_name result threads hashsz additional

ref=$1
ref_base=$2
reads=$3
reads_base=$4
out=$5
threads=$6
hashsz=$7
additional=$8

dir=$(dirname $(readlink -e $0))

echo "build started" 
date 
echo -e "0\tBuilding index\t"`date +%s` >> $out/mapping_time.log
MosaikBuild -fr $ref -oa $out/$ref_base.dat && \
date && \
MosaikBuild -q $reads -out $out/$reads_base.dat -st Illumina && \
echo "align started" && \
date && \
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
MosaikAligner -p $threads -hs $hashsz -in $out/$reads_base.dat -out $out/$reads_base -ia $out/$ref_base.dat -annpe $dir/../../../misc/pe.ann -annse $dir/../../../misc/se.ann && \
echo "align done" && \
date && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log
