#!/bin/bash
#USAGE ./bowtie2.sh reference ref_base_name reads reads_base_name result threads hashsz additional

ref=$1
ref_base=$2
reads=$3
reads_base=$4
out=$5
threads=$6
hashsz=$7
additional=$8

echo "build started" | tee -a $out/log 
date 
echo -e "0\tBuilding index\t"`date +%s` >> $out/mapping_time.log
bowtie2-build $ref $out/$ref_base.ind && \
echo "mapping started" && \
date && \
echo -e "2\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
bowtie2 -t -q -p $threads --sensitive --sam-rg -x $out/$ref_base.ind -U $reads > $out/$reads_base.sam && \
echo "mapping done" && \
date && \
echo -e "3\tAlignment done\t"`date +%s` >> $out/mapping_time.log
