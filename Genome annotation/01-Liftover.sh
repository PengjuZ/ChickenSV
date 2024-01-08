#!/bin/sh
###configuration
#conda create -n GAnno
#conda install -c bioconda liftoff

source /share/apps/anaconda3/bin/activate GAnno
Ref=$1
Refanno=$2
Genome=$3
OUT=$4
liftoff $Genome $Ref -g $Refanno -o $OUT/liftover.gff -u $OUT/unmapped_features.txt -dir $OUT




