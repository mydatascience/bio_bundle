#!/bin/bash
#usage ./mpileup.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

samtools faidx $ref
date
samtools mpileup -uf $ref $bam | bcftools view -cvg - > $out/$name.vcf
date
