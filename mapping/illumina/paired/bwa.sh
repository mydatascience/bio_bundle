#!/bin/bash
#USAGE ./bwa.sh reference ref_base_name reads1 reads1_base_name reads2 reads2_base_name result threads hashsz additional

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

echo "build index started"
date
echo -e "0\tBuilding index\t"`date +%s` >> $out/mapping_time.log
bwa index -a bwtsw $ref && \
echo "align started" && \
date && \
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
bwa aln -t $threads $ref $reads1 > $out/$reads1_base.sai && \
bwa aln -t $threads $ref $reads2 > $out/$reads2_base.sai && \
echo "making sam from sai started" && \
date && \
bwa sampe $ref $out/$reads1_base.sai $out/$reads2_base.sai $reads1 $reads2 > $out/$reads1_base.sam && \
echo "making sam from sai done" && \
date && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log
