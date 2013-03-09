#!/bin/bash
#usage ./mpileup.sh ref.fa aln.bam res_name out_dir bed

ref=$1
bam=$2
name=$3
out=$4
bed=$5

samtools faidx $ref
echo "mpileup started"
date
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log 
if [ -n $bed ]; then
    samtools mpileup -M 30 -q 30 -Q 30 -f $ref -l $bed $bam > $out/$name.bcf
else
    samtools mpileup -M 30 -q 30 -Q 30 -f $ref $bam > $out/$name.bcf
fi
if [ -s $out/$name.bcf ]; then
    echo "SNP search started"
    date
    java -jar $TOOLS_PATH/VarScan.v2.3.4.jar mpileup2indel $out/$name.bcf --min_avg_qual 30 --output-vcf > $out/$name.indel.vcf
    echo "indel search started"
    date
    java -jar $TOOLS_PATH/VarScan.v2.3.4.jar mpileup2snp $out/$name.bcf --min_avg_qual 30 --output-vcf > $out/$name.snp.vcf
    echo "indel search done"
    date
    echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log 
    vcf-concat $out/$name.*.vcf | vcf-sort > $out/$name.vcf
else 
    date
fi
