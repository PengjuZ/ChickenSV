#!/bin/sh
#######Configuration
#conda create -n Abinitio

source /share/apps/anaconda3/bin/activate Abinitio
Genome=$1
OUT=$2
Thread=$3

#######Supervised-augustus
augustus --species=chicken --gff3=on $Genome > $OUT/Augustus.gff
