#!/bin/bash
#usage ./mpileup.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

samtools faidx $ref
echo "mpileup started"
date
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log 
samtools mpileup -f $ref $bam > $out/$name.bcf
if [ -s $out/$name.bcf ]; then
    echo "SNP search started"
    date
    java -jar $TOOLS_PATH/VarScan.v2.3.4.jar mpileup2indel $out/$name.bcf --output-vcf > $out/$name.indel.vcf
    echo "indel search started"
    date
    java -jar $TOOLS_PATH/VarScan.v2.3.4.jar mpileup2snp $out/$name.bcf --output-vcf > $out/$name.snp.vcf
    echo "indel search done"
    date
    echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log 
else 
    date
    vcf-concat $out/$name.*.vcf > $out/$name.vcf
fi
