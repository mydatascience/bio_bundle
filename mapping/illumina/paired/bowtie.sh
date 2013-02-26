#!/bin/bash
#USAGE ./bowtie.sh reference ref_base_name reads1 reads1_base_name reads2 reads2_base_name result threads hashsz additional

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
echo -e "0\tBuilding index\t"`date +%s` >> $out/mapping_time.log
bowtie-build $ref $out/$ref_base.ind && \
echo "mapping started" && \
date && \
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
bowtie -t -q -p $threads -a --fr --minins 0 --maxins 500 --sam $out/$ref_base.ind -1 $reads1 -2 $reads2 > $out/$reads1_base.sam && \
echo "mapping done" && \
date && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log
