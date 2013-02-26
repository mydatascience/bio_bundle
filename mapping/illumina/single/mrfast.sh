#!/bin/bash
#USAGE ./mrfast.sh reference ref_base_name reads reads_base_name result threads hashsz additional

ref=$1
ref_base=$2
reads=$3
reads_base=$4
out=$5
threads=$6
hashsz=$7
additional=$8

echo "build index started" 
date
mrfast --index $ref --ws 12 && \
echo "align started" && \
date && \
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
mrfast --search $ref $additional --seq $reads -o $out/$reads_base.sam && \
echo "align date" && \
date
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log
