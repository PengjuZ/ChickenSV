#!/bin/sh
################################nextDenovo
source /share/apps/anaconda3/bin/activate nextDenovo
realpath $Data/ONT.gz > $OUT/run.fofn
cp ./00-Software/NextDenovo/doc/run.cfg $OUT/run.cfg
./NextDenovo/nextDenovo $OUT/run.cfg

