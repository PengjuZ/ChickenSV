#!/bin/sh

Script1=./00-Software/VGP/Salsa/filter_five_end.pl
Script2=./00-Software/VGP/Salsa/two_read_bam_combiner.pl
OUTA=./03-Cgenome/01-Assembly/04-HIC/01-Pacbio/Hap1
IN1=./03-Cgenome/00-Data/HIC_1.QC.fq.gz
IN2=./03-Cgenome/00-Data/HIC_2.QC.fq.gz
Data=./03-Cgenome/00-Data

mkdir $OUTA/Salsa2
OUT=$OUTA/Salsa2
source /share/apps/anaconda3/bin/activate Ntools
cat ./03-Cgenome/01-Assembly/01-Pacbio-trioHIC/Pacbio.hic.hap1.fa > $OUT/Chicken.fa
hap=$OUT/Chicken.fa
bwa index $hap
samtools faidx $hap
bwa mem -t 28 $hap $IN1 | samtools view -Sb - > $OUT/hap.HIC_1.bam
samtools view -h $OUT/hap.HIC_1.bam | perl $Script1 | samtools view -@28 -Sb - > $OUT/hap.F.HIC_1.bam
bwa mem -t 28 $hap $IN2 | samtools view -Sb - > $OUT/hap.HIC_2.bam
samtools view -h $OUT/hap.HIC_2.bam | perl $Script1 | samtools view -@28 -Sb - > $OUT/hap.F.HIC_2.bam
perl $Script2 $OUT/hap.F.HIC_1.bam $OUT/hap.F.HIC_2.bam | samtools view -@28 -Sb > $OUT/hap.HIC.bam
samtools sort -@28 -m2G -O bam -o $OUT/hap.HIC.sort.bam $OUT/hap.HIC.bam
samtools index $OUT/hap.HIC.sort.bam
picard -Xmx280g MarkDuplicates --INPUT $OUT/hap.HIC.sort.bam --OUTPUT $OUT/hap.bam --METRICS_FILE $OUT/hap.metrics.txt --CREATE_INDEX true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 --REMOVE_DUPLICATES true --TMP_DIR $OUT/00-TMP/ --VALIDATION_STRINGENCY LENIENT --SORTING_COLLECTION_SIZE_RATIO 0.1 --MAX_RECORDS_IN_RAM 500000
samtools sort -@28 -n -m2800m -O BAM -o $OUT/hap.dedup.bam $OUT/hap.bam
bedtools bamtobed -i $OUT/hap.dedup.bam > $OUT/hap.dedup.bed

source /share/apps/anaconda3/bin/activate Salsa2
./00-Software/SALSA-2.3/run_pipeline.py -a $hap -l $hap.fai -e GATC -b $OUT/hap.dedup.bed -o $OUT/ -m yes -p yes 


mkdir $OUTA/YaHS
OUT=$OUTA/YaHS
source /share/apps/anaconda3/bin/activate Ntools
cat $OUTA/Salsa2/assembly.cleaned.fasta > $OUT/Chicken.fa
hap=$OUT/Chicken.fa
bwa index $hap
samtools faidx $hap
bwa mem -t 28 $hap $IN1 | samtools view -Sb - > $OUT/hap.HIC_1.bam
samtools view -h $OUT/hap.HIC_1.bam | perl $Script1 | samtools view -@28 -Sb - > $OUT/hap.F.HIC_1.bam
bwa mem -t 28 $hap $IN2 | samtools view -Sb - > $OUT/hap.HIC_2.bam
samtools view -h $OUT/hap.HIC_2.bam | perl $Script1 | samtools view -@28 -Sb - > $OUT/hap.F.HIC_2.bam
perl $Script2 $OUT/hap.F.HIC_1.bam $OUT/hap.F.HIC_2.bam | samtools view -@28 -Sb > $OUT/hap.HIC.bam
samtools sort -@28 -m2G -O bam -o $OUT/hap.HIC.sort.bam $OUT/hap.HIC.bam
samtools index $OUT/hap.HIC.sort.bam
picard -Xmx280g MarkDuplicates --INPUT $OUT/hap.HIC.sort.bam --OUTPUT $OUT/hap.bam --METRICS_FILE $OUT/hap.metrics.txt --CREATE_INDEX true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 --REMOVE_DUPLICATES true --TMP_DIR $OUT/00-TMP/ --VALIDATION_STRINGENCY LENIENT --SORTING_COLLECTION_SIZE_RATIO 0.1 --MAX_RECORDS_IN_RAM 500000
samtools sort -@28 -n -m2800m -O BAM -o $OUT/hap.dedup.bam $OUT/hap.bam

source /share/apps/anaconda3/bin/activate YaHS
./00-Software/yahs-1.2a/yahs -o $OUT/ChickenHIC -e GATC -q 10 $hap $OUT/hap.dedup.bam



