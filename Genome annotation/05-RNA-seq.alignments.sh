#!/bin/sh
#conda create -n Braker2
#conda install -c anaconda perl
#conda install -c bioconda perl-app-cpanminus
#cpanm YAML Hash::Merge Logger::Simple Parallel::ForkManager
#conda install -c bioconda perl-scalar-util-numeric perl-class-data-inheritable
#conda install -c bioconda perl-exception-class perl-test-pod perl-file-which perl-mce perl-threaded perl-list-util perl-math-utils 
#conda install -c anaconda biopython
#conda install -c bioconda cdbtools
#cpanm File::HomeDir List::MoreUtils Hash::Merge
#conda install -c bioconda perl-test-leaktrace
#conda install -c bioconda augustus=3.4.0
#conda install -c anaconda python=3.8.0
#conda install -c anaconda samtools=1.7
#conda install -c bioconda spaln
#conda install -c bioconda exonerate
#conda install -c bioconda vardict-java
#conda install -c bioconda hisat2

GenomeM=$1
GenomeUM=$2
OUTdir=$3
RNAlist=$4
Thread=$5
Data=./03-Cgenome/00-Data/Raw/00-Ready/03-RNA-seq/

cp ./00-Software/GeneMark/gm_key_64 ~/.gm_key
source /share/apps/anaconda3/bin/activate Braker2
export GENEMARK_PATH=./00-Software/GeneMark/gmes_linux_64_4
export GUSHR_PATH=./00-Software/GeneMark/GUSHR-master
PERL5LIB=$PERL5LIB:./.conda/envs/GeneMark/lib/site_perl/5.26.2; export PERL5LIB

:<<Masknote
Masknote
mkdir $OUTdir/Index
cd $OUTdir/Index
cat $GenomeM > GenomeM.fa
hisat2-build GenomeM.fa $OUTdir/Index/GenomeM

mkdir $OUTdir/Mapping
cd $OUTdir/Mapping

for Sample in $(cat $RNAlist | sed 's/\t/:/g')
do
SRA=(${Sample//\:/ })
./00-Software/bbmap/bbduk.sh -Xmx16g in1=$Data/${SRA[0]}\_1.QC.fq.gz in2=$Data/${SRA[0]}\_2.QC.fq.gz out1=$OUTdir/Mapping/${SRA[0]}\_1.fq.gz out2=$OUTdir/Mapping/${SRA[0]}\_2.fq.gz minlen=70 qtrim=rl trimq=10 ktrim=r k=23 mink=11 ref=./00-Software/bbmap/resources/truseq.fa.gz hdist=1
hisat2 -p $Thread -x $OUTdir/Index/GenomeM -1 $OUTdir/Mapping/${SRA[0]}\_1.fq.gz -2 $OUTdir/Mapping/${SRA[0]}\_2.fq.gz | samtools sort -@ $Thread > $OUTdir/Mapping/${SRA[0]}.sorted.bam
samtools index $OUTdir/Mapping/${SRA[0]}.sorted.bam
done

samtools merge -@ $Thread $OUTdir/Merge.bam `ls *.sorted.bam`

mkdir $OUTdir/OUT

perl ./00-Software/BRAKER-2.1.6/scripts/braker.pl --genome=$GenomeM --species=chicken --bam=$OUTdir/Merge.bam --workingdir=$OUTdir/OUT --cores $Thread --UTR=on --gff3 --softmasking --useexisting
#


