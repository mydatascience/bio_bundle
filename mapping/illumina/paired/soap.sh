#!/bin/bash
#USAGE ./bowtie2.sh reference ref_base_name reads1 reads1_base_name reads2 reads2_base_name result threads hashsz additional

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
2bwt-builder $ref && \
echo "mapping started" && \
date && \
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
soap -a $reads1 -2 $reads2 -D $ref.index -m 0 -x 500 -o $out/$reads1_base.soap -p $threads -g 1 && \
echo "mapping done" && \
date && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log
$dir/../../../misc/soap2sam.pl $out/$reads1_base.soap > $out/$reads1_base.sam
