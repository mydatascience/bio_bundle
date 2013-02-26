#!/bin/bash
#usage ./gatk.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

java -jar $TOOLS_PATH/picard-tools-1.84/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=SILENT RGLB=$ref RGPL=Unknown RGSM=unknown INPUT=$bam OUTPUT=$out/$name.rg.bam RGPU=run && \
samtools index $out/$name.rg.bam && \
date && \
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log && \
java -jar $TOOLS_PATH/GenomeAnalysisTKLite.jar -R $ref -I $out/$name.rg.bam -T UnifiedGenotyper -o $out/$name.vcf && \
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log && \
date && \
rm $out/$name.rg.bam
