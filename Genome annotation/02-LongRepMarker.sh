#!/bin/sh
source /share/apps/anaconda3/bin/activate TEanno212
Genome=$1
OUT=$2
cp -r ./00-Software/LongRepMarker_v2.1.2-master/ $OUT/LongRepMarker/
cd $OUT/LongRepMarker/
java -Xmx1000g LongRepMarker -r $Genome -T no -k 49 -m 100 -t 10 -R fast -o $OUT/
