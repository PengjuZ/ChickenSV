#!/bin/sh

Hap1=./Pacbio.hic.hap1.fa
Hap2=./Pacbio.hic.hap2.fa
source /share/apps/anaconda3/bin/activate Ntools
mkdir $WorkDir/Split
cd $WorkDir/Split
bwa index $Hap1
bwa index $Hap2
bwa mem $Hap1 $Data/NGS_1.QC.fq.gz $Data/NGS_2.QC.fq.gz \
    | samtools view -bh - | samtools sort -n - > Hap1_ref.bam
bwa mem $Hap2 $Data/NGS_1.QC.fq.gz $Data/NGS_2.QC.fq.gz \
    | samtools view -bh - | samtools sort -n - > Hap2_ref.bam
source /share/apps/anaconda3/bin/activate verkko
$Software/classify_by_alignment --hapA-in Hap1_ref.bam \
    --hapA-out Hap1_classified.bam \
    --hapB-in Hap2_ref.bam \
    --hapB-out Hap2_classified.bam
#bam2fastq
source /share/apps/anaconda3/bin/activate Ntools
bwa mem $Hap1 $Data/HIC_1.QC.fq.gz $Data/HIC_2.QC.fq.gz \
    | samtools view -bh - | samtools sort -n - > Hap1.HIC_ref.bam
bwa mem $Hap2 $Data/HIC_1.QC.fq.gz $Data/HIC_2.QC.fq.gz \
    | samtools view -bh - | samtools sort -n - > Hap2.HIC_ref.bam
conda /share/apps/anaconda3/bin/activate verkko
$Software/classify_by_alignment --hapA-in Hap1.HIC_ref.bam \
    --hapA-out Hap1.HIC_classified.bam \
    --hapB-in Hap2.HIC_ref.bam \
    --hapB-out Hap2.HIC_classified.bam

source /share/apps/anaconda3/bin/activate Ntools
minimap2 -H -d ./Pacbio.hic.hap1.mmi ./Pacbio.hic.hap1.fa
minimap2 -H -d ./Pacbio.hic.hap2.mmi ./Pacbio.hic.hap2.fa
minimap2 -ax map-hifi $Hap1 $Data/HIFI.filt.fastq.gz | samtools view -bh - | samtools sort -n - > Hap1.HIFI_ref.bam
minimap2 -ax map-hifi $Hap2 $Data/HIFI.filt.fastq.gz | samtools view -bh - | samtools sort -n - > Hap2.HIFI_ref.bam
#conda /share/apps/anaconda3/bin/activate verkko
$Software/classify_by_alignment --hapA-in Hap1.HIFI_ref.bam \
    --hapA-out Hap1.HIFI_classified.bam \
    --hapB-in Hap2.HIFI_ref.bam \
    --hapB-out Hap2.HIFI_classified.bam

source /share/apps/anaconda3/bin/activate Ntools
minimap2 -ax map-ont $Hap1 /share/home/zju_zhaopj/03-Cgenome/03-ONT-Ratatosk/ONT.correct.fastq | samtools view -bh - | samtools sort -n - > Hap1.ONT_ref.bam
minimap2 -ax map-ont $Hap2 /share/home/zju_zhaopj/03-Cgenome/03-ONT-Ratatosk/ONT.correct.fastq | samtools view -bh - | samtools sort -n - > Hap2.ONT_ref.bam
conda /share/apps/anaconda3/bin/activate verkko
$Software/classify_by_alignment --hapA-in Hap1.ONT_ref.bam \
    --hapA-out Hap1.ONT_classified.bam \
    --hapB-in Hap2.ONT_ref.bam \
    --hapB-out Hap2.ONT_classified.bam


