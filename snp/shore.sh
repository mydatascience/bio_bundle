#!/bin/bash
#usage ./shore.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

echo "preprocess started"
date
shore preprocess -f $ref -i $out/index && \
shore preprocess --upgrade -g -U $out/index/*.fasta.shore && \
echo "converting bam to maplist" && \
date && \
shore convert -r $ref Alignment2Maplist $bam > $out/map.list && \
echo "start snp calling" && \
date && \
shore qVar -n Unknown -f $out/index/*.fasta.shore -i $out/map.list -s $SHORE_SCORE/scoring_matrix_het.txt -o $out/res && \
echo "snp calling done" && \
date
