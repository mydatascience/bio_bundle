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

dir=$(dirname $(readlink -e $0))

echo "build started" 
date 
echo -e "0\tBuilding index\t"`date +%s` >> $out/mapping_time.log
2bwt-builder $ref && \
echo "mapping started" && \
date && \
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
soap -a $reads -D $ref.index -o $out/$reads_base.soap -p $threads -g 1 && \
echo "mapping done" && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log && \
date 
$dir/../../../misc/soap2sam.pl $out/$reads_base.soap > $out/$reads_base.sam
