#!/bin/bash
#usage ./gatk.sh ref.fa aln.bam res_name out_dir bed

ref=$1
bam=$2
name=$3
out=$4
bed=$5

java -jar $TOOLS_PATH/picard-tools-1.84/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=SILENT RGLB=$ref RGPL=Unknown RGSM=unknown INPUT=$bam OUTPUT=$out/$name.rg.bam RGPU=run && \
samtools index $out/$name.rg.bam && \
date && \
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log && \
if [ -n $bed ]; then
#    java -jar $TOOLS_PATH/GenomeAnalysisTKLite.jar -R $ref -L $bed -U ALLOW_SEQ_DICT_INCOMPATIBILITY -I $out/$name.rg.bam -T UnifiedGenotyper -o $out/$name.vcf
#    java -jar $TOOLS_PATH/GenomeAnalysisTKLite.jar -R $ref -L $bed -mbq 20 -dcov 2000 -deletions 2.0 -stand_call_conf 30 -stand_emit_conf 10 -glm INDEL -U ALLOW_SEQ_DICT_INCOMPATIBILITY -I $out/$name.rg.bam -T UnifiedGenotyper -o $out/$name.vcf
    java -jar $TOOLS_PATH/GenomeAnalysisTKLite.jar -minIndelCnt 1 -R $ref -L $bed -stand_call_conf 30 -stand_emit_conf 30 -glm INDEL -A DepthOfCoverage \
        -A HaplotypeScore \
            -A InbreedingCoeff \
            -U ALLOW_SEQ_DICT_INCOMPATIBILITY -I $out/$name.rg.bam -T UnifiedGenotyper -o $out/$name.vcf
else
    java -jar $TOOLS_PATH/GenomeAnalysisTKLite.jar -R $ref -I $out/$name.rg.bam -T UnifiedGenotyper -o $out/$name.vcf
fi
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log
date
rm $out/$name.rg.bam
