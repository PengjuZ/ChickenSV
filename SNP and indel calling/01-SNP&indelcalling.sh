#!/bin/bash
####################################################################################
###Files
ID=$1
Index=$2
Outdir=$3
Thread=$4

Datadir=./03-Cgenome/00-Data/Population/Clean_data
###SoftWares
module load anaconda3/4.12.0

fastp=./00-Software/fastp

###
mkdir $Outdir/$ID/
mkdir $Outdir/$ID/Run
RUN=$Outdir/$ID/Run
mkdir $Outdir/$ID/Out
OUT=$Outdir/$ID/Out

cat $Datadir/$ID\_clean_1.fq.gz > $RUN/Sample_1.fastq.gz
cat $Datadir/$ID\_clean_2.fq.gz > $RUN/Sample_2.fastq.gz
$fastp -i $RUN/Sample_1.fastq.gz -I $RUN/Sample_2.fastq.gz -o $RUN/Sample_1.QC.fastq.gz -O $RUN/Sample_2.QC.fastq.gz -j $OUT/$ID.fastp.json -h $OUT/$ID.fastp.html -q 20 -u 30 -l 75 -w $Thread
Fastq1=$RUN/Sample_1.QC.fastq.gz
Fastq2=$RUN/Sample_2.QC.fastq.gz

bwa mem -t $Thread -R "@RG\tID:$ID\tLB:$ID\tPL:ILLUMINA\tSM:$ID" -M $Index/Genome.fa $Fastq1 $Fastq2 | sentieon util sort -r $Index/Genome.fa -o $RUN/$ID.bam -t $Thread --sam2bam -i -

sentieon driver -t $Thread --temp_dir $Outdir/TMP/ -i $RUN/$ID.bam --algo LocusCollector --fun score_info $RUN/$ID.score.txt
sentieon driver -t $Thread --temp_dir $Outdir/TMP/ -i $RUN/$ID.bam --algo Dedup --rmdup --score_info $RUN/$ID.score.txt --metrics $RUN/$ID.rmdup_metrics.txt $RUN/$ID.rmdup.bam

sentieon driver -r $Index/Genome.fa -t $Thread --temp_dir $Outdir/TMP/ -i $RUN/$ID.rmdup.bam --algo Realigner $OUT/$ID.realigner.bam
sentieon driver -r $Index/Genome.fa -t $Thread --temp_dir $Outdir/TMP/ -i $OUT/$ID.realigner.bam --algo QualCal $RUN/$ID.data.table

sentieon driver -r $Index/Genome.fa -t $Thread --temp_dir $Outdir/TMP/ -i $OUT/$ID.realigner.bam -q $RUN/$ID.data.table --algo Haplotyper --emit_conf=30 --call_conf=30 --emit_mode gvcf --genotype_model multinomial $OUT/$ID.g.vcf.gz

samtools stats $OUT/$ID.realigner.bam -@ $Thread > $OUT/$ID.stats.txt
samtools flagstat $OUT/$ID.realigner.bam -@ $Thread > $OUT/$ID.flagstat.txt


