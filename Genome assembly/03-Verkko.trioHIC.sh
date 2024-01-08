#!/bin/sh
WorkDir=./02-Hybrid-trioHIC
Threads="40"

Hap1=./01-Pacbio-trioHIC/Pacbio.hic.hap1.fa
Hap2=./01-Pacbio-trioHIC/Pacbio.hic.hap2.fa

source /share/apps/anaconda3/bin/activate Verkko
mkdir $WorkDir/meryl/
./verkko-v1.1/lib/verkko/bin/meryl count k=21 threads=$Threads compress $Hap1 output $WorkDir/meryl/hap1.meryl
./verkko-v1.1/lib/verkko/bin/meryl count k=21 threads=$Threads compress $Hap2 output $WorkDir/meryl/hap2.meryl
./verkko-v1.1/lib/verkko/bin/meryl difference k=21 threads=$Threads $WorkDir/meryl/hap1.meryl $WorkDir/meryl/hap2.meryl output $WorkDir/meryl/hap1.meryl.only
./verkko-v1.1/lib/verkko/bin/meryl difference k=21 threads=$Threads $WorkDir/meryl/hap2.meryl $WorkDir/meryl/hap1.meryl output $WorkDir/meryl/hap2.meryl.only
./verkko-v1.1/bin/verkko -d $WorkDir --hifi $Data/HIFI.filt.fastq.gz --nano $Data/ONT.fastq.gz --threads $Threads --hap-kmers $WorkDir/meryl/hap1.meryl.only $WorkDir/meryl/hap2.meryl.only hic
