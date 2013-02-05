#!/bin/bash
#usage ./gatk.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

picard-tools AddOrReplaceReadGroups RGLB=$ref RGPL=Unknown RGSM=unknown INPUT=$reads OUTPUT=$out/aln_rg.bam RGPU=run &&
date
GenomeAnalysisTKLite.jar -R $reads -I $out/aln_rg.bam -T UnifiedGenotyper -o $out/out.vcf
date
