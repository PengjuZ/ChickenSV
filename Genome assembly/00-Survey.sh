#!/bin/bash
source /share/apps/anaconda3/bin/activate genomescope2

fq1=./01-Assembly/00-Survey/NGS_1.QC.fq
fq2=./01-Assembly/00-Survey/NGS_2.QC.fq
jellyfish count -C -m 21 -s 1000000000 -o ./01-Assembly/00-Survey/NGS.jf -t 15 ${fq1} ${fq2}
jellyfish histo -t 15 -o ./01-Assembly/00-Survey/NGS.histo ./01-Assembly/00-Survey/NGS.jf
genomescope2 -i ./01-Assembly/00-Survey/NGS.histo -o ./01-Assembly/00-Survey/NGS.Result -k 21 -p 2
