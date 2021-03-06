#!/bin/bash
#USAGE ./stampy.sh reference ref_base_name reads reads_base_name result threads hashsz additional

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
echo -e "0\tBuilding index\t"`date +%s` >> $out/mapping_time.log
python $StampyPATH/stampy.py --noparseNCBI -G $out/$ref_base $ref && \
echo "building hash started" && \
date && \
python $StampyPATH/stampy.py -g $out/$ref_base -H $out/$ref_base && \
echo "align started" && \
date && \
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
python $StampyPATH/stampy.py $additional -g $out/$ref_base -h $out/$ref_base -t $threads --bamkeepgoodreads -M $reads > $out/$reads_base.sam && \
echo "align done" && \
date && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log
