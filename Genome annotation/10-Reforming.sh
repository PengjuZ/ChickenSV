#!/bin/sh

######
GenomeUM=$1
Workdir=$2
ID=$3
Thread=$4
liftovergff=$5
PASAgff=$(ls $Workdir/09-Asplicing/*.gff3 | grep "updates")
######
source /share/apps/anaconda3/bin/activate PopVa
:<<Masknote
Masknote
##01-data
cat $GenomeUM > 00-$ID.fa
cat $liftovergff > 01-liftover.gff
cat $PASAgff > 01-PASAgff

##02-QC
singularity run agat_1.0.0--pl5321hdfd78af_0.sif
agat config --expose
agat_convert_sp_gxf2gxf.pl --gff 01-liftover.gff -o 02-liftover.gff
agat_convert_sp_gxf2gxf.pl --gff 01-PASAgff -o 02-PASAgff

##03-Merge
agat_sp_complement_annotations.pl --ref 02-liftover.gff --add 02-PASAgff --out 03-Merge.gff

##04-QC&Filter
source /share/apps/anaconda3/bin/activate PopVa
gff3_QC -g 03-Merge.gff -f 00-$ID.fa -o 04-sample.qc -s 04-statistic.txt
gff3_fix -qc_r 04-sample.qc -g 03-Merge.gff -og 04-Merge.gff
