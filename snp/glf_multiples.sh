#!/bin/bash
#usage ./glf_multiples.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

echo "pileup started"
date
samtools-hybrid pileup -gf $ref $bam > $out/$name.glf &&
echo "snp calling started"
date
glfMultiples -b $out/$name.vcf $out/$name.glf
echo "snp calling done"
date
