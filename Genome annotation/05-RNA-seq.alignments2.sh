#!/bin/sh
#conda create -n RNAssembly
#conda install -c bioconda trinity
#conda install -c anaconda transdecoder
#conda install -c bioconda hisat2
#conda install -c bioconda cufflinks

GenomeM=$1
GenomeUM=$2
OUTdir=$3
RNAlist=$4
Thread=$5
Datas=$6
Data=$Datas/Mapping

source /share/apps/anaconda3/bin/activate RNAssembly

mkdir $OUTdir/Index
cd $OUTdir/Index
cat $GenomeUM > GenomeUM.fa
hisat2-build GenomeUM.fa $OUTdir/Index/GenomeUM
mkdir $OUTdir/Mapping
cd $OUTdir/Mapping

for Sample in $(cat $RNAlist | sed 's/\t/:/g')
do
SRA=(${Sample//\:/ })
hisat2 --dta -p $Thread -x $OUTdir/Index/GenomeUM -1 $Data/${SRA[0]}\_1.fq.gz -2 $Data/${SRA[0]}\_2.fq.gz | samtools sort -@ $Thread > $OUTdir/Mapping/${SRA[0]}.sorted.bam
samtools index $OUTdir/Mapping/${SRA[0]}.sorted.bam
done

samtools merge -@ $Thread -o $OUTdir/Merge.bam `ls *.sorted.bam`
stringtie -p $Thread -o $OUTdir/Merge.gtf $OUTdir/Merge.bam

mkdir $OUTdir/OUT
cd $OUTdir/OUT

gtf_genome_to_cdna_fasta.pl $OUTdir/Merge.gtf $OUTdir/Index/GenomeUM.fa > transcripts.fasta
gtf_to_alignment_gff3.pl $OUTdir/Merge.gtf > transcripts.gff3
TransDecoder.LongOrfs -t transcripts.fasta
TransDecoder.Predict -t transcripts.fasta
cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     transcripts.gff3 \
     transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
