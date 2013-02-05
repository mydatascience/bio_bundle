#!/bin/bash
#USAGE ./bowtie.sh reference ref_base_name reads reads_base_name result threads hashsz additional

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
bwa index -a bwtsw $ref && \
echo "align started" && \
date && \
bwa aln -t $threads $ref $reads > $out/$reads_base.sai && \
echo "making sam from sai started" && \
date && \
bwa samse $ref $out/$reads_base.sai $reads > $out/$reads_base.sam && \
echo "making sam from sai done" && \
date
