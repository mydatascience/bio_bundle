#!/bin/bash
#usage ./shore.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

dir=$(dirname $(readlink -e $0))

echo "preprocess started"
date
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log 
shore preprocess -f $ref -i $out/index && \
shore preprocess --upgrade -g -U $out/index/$(basename $ref).shore && \
echo "converting bam to maplist" && \
date && \
shore convert -r $ref Alignment2Maplist $bam > $out/map.list && \
echo "start snp calling" && \
date && \
shore qVar -n Unknown -f $out/index/$(basename $ref).shore -i $out/map.list -s $dir/../misc/scoring_matrix_het.txt -o $out/$name && \
echo "snp calling done" && \
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log && \
date && \
vcf-concat $out/$name/ConsensusAnalysis/*.vcf > $out/$name.vcf && \
echo "succeed" && \
date
