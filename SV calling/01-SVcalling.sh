#!/bin/sh
###configuration
#conda create -n SVcalling
#conda install -c bioconda manta
#conda install -c bioconda delly
#conda install -c bioconda gridss
#conda install -c bioconda wham
#conda install -c bioconda smoove
#conda create -n Dysgu python=3.7
#pip install numpy dysgu


ID=$1
Index=$2
Outdir=$3
Thread=$4

Ref=$Index/Genome.fa
Bam=$Outdir/$ID/Out/$ID.realigner.bam
mkdir $Outdir/$ID/Out/SVOUT
SVOUT=$Outdir/$ID/Out/SVOUT

###Manta (RP+SR+AS)
source /share/apps/anaconda3/bin/activate SVcalling
mkdir $Outdir/$ID/Run/Manta
cd $Outdir/$ID/Run/Manta
configManta.py --bam $Bam --referenceFasta $Ref --runDir $Outdir/$ID/Run/Manta
python $Outdir/$ID/Run/Manta/runWorkflow.py
gunzip $Outdir/$ID/Run/Manta/results/variants/candidateSV.vcf.gz
cat $Outdir/$ID/Run/Manta/results/variants/candidateSV.vcf > $SVOUT/Manta.vcf

###Delly (RP+SR)
mkdir $Outdir/$ID/Run/Delly
cd $Outdir/$ID/Run/Delly
delly call -o delly.bcf -g $Ref $Bam
bcftools view delly.bcf > delly.vcf
cat $Outdir/$ID/Run/Delly/delly.vcf > $SVOUT/Delly.vcf

###Wham (RP+SR)
mkdir $Outdir/$ID/Run/Wham
cd $Outdir/$ID/Run/Wham
whamg -x $Thread -a $Ref -f $Bam > Wham.vcf 2> Wham.err
cat $Outdir/$ID/Run/Wham/Wham.vcf > $SVOUT/Wham.vcf

###Smoove (lumpy:RP+SR+RD)
mkdir $Outdir/$ID/Run/Smoove
cd $Outdir/$ID/Run/Smoove
smoove call --outdir $Outdir/$ID/Run/Smoove/ --name $ID --fasta $Ref -p $Thread --genotype $Bam
gunzip $Outdir/$ID/Run/Smoove/$ID\-smoove.genotyped.vcf.gz
cat $Outdir/$ID/Run/Smoove/$ID\-smoove.genotyped.vcf > $SVOUT/Smoove.vcf

###Dysgu
mkdir $Outdir/$ID/Run/Dysgu
cd $Outdir/$ID/Run/Dysgu
source /share/apps/anaconda3/bin/activate Dysgu
dysgu run -p $Thread $Ref $Outdir/$ID/Run/Dysgu/tmp $Bam > Dysgu.vcf
cat $Outdir/$ID/Run/Dysgu/Dysgu.vcf > $SVOUT/Dysgu.vcf

###GRIDSS2 (RP+SR+AS)
source /share/apps/anaconda3/bin/activate SVcalling
export PATH=./00-Software/samtools-1.16.1/bin:$PATH
mkdir $Outdir/$ID/Run/Gridss2
cd $Outdir/$ID/Run/Gridss2
./00-Software/Gridss/gridss --reference $Ref --output Gridss2.vcf --threads $Thread --jar ./00-Software/Gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar --workingdir $Outdir/$ID/Run/Gridss2 $Bam
cat $Outdir/$ID/Run/Gridss2/Gridss2.vcf > $SVOUT/Gridss2.vcf

###Survivor
cd $Outdir/$ID/Run/
samtools view -H $Bam > lowMQ.sam
samtools view $Bam | awk '$5<5 {print $0}' >> lowMQ.sam
samtools view -S -b -h lowMQ.sam > lowMQ.bam
samtools depth lowMQ.bam > lowMQ.cov
SURVIVOR bincov lowMQ.cov 10 2 > $Outdir/$ID/Out/$ID.lowMQ.bed

cd $SVOUT/
ls $SVOUT/*vcf > Sample.txt
SURVIVOR merge Sample.txt 1000 2 1 1 0 50 $Outdir/$ID/Out/$ID.SV.vcf


