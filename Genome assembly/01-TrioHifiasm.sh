#!/bin/sh

source /share/apps/anaconda3/bin/activate hifiasm
hifiasm -o $OUT/Pacbio -t 24 --h1 $Data/HIC_1.QC.fq.gz --h2 $Data/HIC_2.QC.fq.gz $Data/HIFI.filt.fastq.gz

