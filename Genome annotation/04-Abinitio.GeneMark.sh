#!/bin/sh
#######Configuration
#conda create -n GeneMark
#conda install -c anaconda perl
#conda install -c conda-forge perl-app-cpanminus
#cpanm YAML Hash::Merge Logger::Simple Parallel::ForkManager
#cpanm MCE::Mutex

source /share/apps/anaconda3/bin/activate GeneMark
Genome=$1
OUT=$2
Thread=$3
gmes_petap=/share/home/zju_zhaopj/00-Software/GeneMark/gmes_linux_64_4/gmes_petap.pl
#######Unsupervised-GeneMark-ES
cd $OUT
perl $gmes_petap --ES --cores $Thread --sequence $Genome
gffread genemark.gtf -o- > Genemark.gff
