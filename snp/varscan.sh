#!/bin/bash
#usage ./mpileup.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

samtools faidx $ref
echo "SNP search started"
date
samtools mpileup -f $ref $bam | VarScan.v2.3.4.jar mpileup2indel > $out/$name.indel.vcf
echo "indel search started"
date
samtools mpileup -f $ref $bam | VarScan.v2.3.4.jar mpileup2snp > $out/$name.indel.vcf
date
