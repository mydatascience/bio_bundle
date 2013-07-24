#!/bin/bash
#usage ./gatk.sh ref.fa aln.bam res_name out_dir bed

ref=$1
bam=$2
name=$3
out=$4
bed=$5

touch $bam.rg
cat $bam | while read line 
do
    basename=$(basename $line)
    java -jar $TOOLS_PATH/picard-tools-1.84/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=SILENT RGLB=$ref RGPL=Unknown RGSM=$line INPUT=$line OUTPUT=$out/$basename.rg.bam RGPU=run && \
    samtools index $out/$basename.rg.bam && \
    echo "$out/$basename.rg.bam" >> $bam.rg
done
date && \
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log && \
if [ -n "$bed" ]; then
    cat $bam.rg | sed 'N;s/\n/ -I /;' | xargs -I {} -n 1 -P 1 sh -c "java -jar $TOOLS_PATH/GenomeAnalysisTK.jar -R $ref -L $bed -dcov 1000 -deletions 1 -mbq 0 -stand_call_conf 0 -stand_emit_conf 0 -glm BOTH -I {} -T UnifiedGenotyper -o $out/$name.vcf"
else
    cat $bam.rg | sed 'N;s/\n/ -I /;' | xargs -I {} -n 1 -P 1 sh -c "java -jar $TOOLS_PATH/GenomeAnalysisTK.jar -R $ref -dcov 1000 -deletions 1 -mbq 0 -stand_call_conf 0 -stand_emit_conf 0 -glm BOTH -I {} -T UnifiedGenotyper -o $out/$name.vcf"
fi
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log
date
cat $bam.rg | xargs -I {} -n 1 -P 1 sh -c "rm {}"
