#!/bin/sh
#conda create -n PASA
#conda install -c bioconda pasa
#pasa-2.5.3
Pasa_gff3_validator=./.conda/envs/PASA/opt/pasa-2.5.3/misc_utilities/pasa_gff3_validator.pl
Load_Current_Gene_Annotations_dbi=./.conda/envs/PASA/opt/pasa-2.5.3/scripts/Load_Current_Gene_Annotations.dbi

GenomeUM=$1
Workdir=$2
Rundir=$3
Thread=$4
Refanno=$5
RNAdb=$6
ID=$7
Scripts=$8

export PATH=./00-Software/samtools-1.16.1/bin:$PATH
source /share/apps/anaconda3/bin/activate PASA

Launch_PASA_pipeline.pl \
-c $Scripts/09-$ID.align.conf \
-C -R --ALIGNER gmap --CPU $Thread \
-g $GenomeUM \
-t $Workdir/05-RNA-seq/Assembly/OUT/transcripts.fasta

cp $Workdir/08-Merging/Partition/EVM.all.gff EVM.all.gff3
$Pasa_gff3_validator EVM.all.gff3

$Load_Current_Gene_Annotations_dbi \
-c $Scripts/09-$ID.align.conf \
-g $GenomeUM \
-P EVM.all.gff3

Launch_PASA_pipeline.pl \
-c $Scripts/09-$ID.annotCompare.conf \
-A \
--ALIGNER gmap --CPU $Thread \
-g $GenomeUM \
-t $Workdir/05-RNA-seq/Assembly/OUT/transcripts.fasta
