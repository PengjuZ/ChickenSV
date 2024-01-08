#!/bin/sh
#source /share/apps/anaconda3/bin/activate GAnno
#conda create -n RepeatMasker
#conda install -c bioconda repeatmasker
source /share/apps/anaconda3/bin/activate RepeatMasker

LIB=$1
Thread=$2
OUT=$3
Genome=$4
RepeatMasker -lib $LIB -pa $Thread -nolow -gff -dir $OUT/ $Genome
