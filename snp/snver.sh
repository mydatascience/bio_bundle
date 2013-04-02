#!/bin/bash
#usage freebayes.sh ref.fa aln.bam res_name out_dir bed

ref=$1
bam=$2
name=$3
out=$4
bed=$5

date
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log 
if [ -n "$bed" ]; then
    java -jar $TOOLS_PATH/SNVerIndividual.jar -b 0.2 -mq 30 -bq 20 -a 2 -i $bam -r $ref -l $bed -o $out/$name 
else
    java -jar $TOOLS_PATH/SNVerIndividual.jar -b 0.2 -mq 30 -bq 20 -a 2 -i $bam -r $ref -o $out/$name
fi
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log && \
date && \
vcf-concat $out/$name.filter.vcf $out/$name.indel.filter.vcf | vcf-sort > $out/$name.vcf
