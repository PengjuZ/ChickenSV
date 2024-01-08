#!/bin/sh
#conda create -n Homology
#conda install -c bioconda exonerate
#conda install -c bioconda genomethreader
#conda install -c anaconda perl
#conda install -c conda-forge perl-app-cpanminus
#conda install -c bioconda blast
#conda install -c anaconda cd-hit=4.8.1
#conda install -c bioconda bedtools
#conda install seqtk


GenomeM=$1
OUTdir=$2
Thread=$3
ProteinDB1=./03-Cgenome/02-Annotation/00-Proteins/NCBI.fa
ProteinDB2=./03-Cgenome/02-Annotation/00-Proteins/UniProtKB.fa
Exonerate=./00-Software/exonerate-2.3.0-x86_64/bin/exonerate
###CD-hit=./00-Software/cd-hit-v4.6.8-2017-1208/cd-hit

###Protein database
mkdir $OUTdir/ProteinDB
cat $ProteinDB1 $ProteinDB2 > $OUTdir/ProteinDB/ProteinDB.fa
###Remove redundancy
source /share/apps/anaconda3/bin/activate PanCattle
seqkit rmdup -s $OUTdir/ProteinDB/ProteinDB.fa > $OUTdir/ProteinDB/ProteinDB.R.fa
seqtk rename $OUTdir/ProteinDB/ProteinDB.R.fa CHI > $OUTdir/ProteinDB/ProteinDB.RM.fa

###Genome Index
source /share/apps/anaconda3/bin/activate Homology
mkdir $OUTdir/GenomeIndex
cd $OUTdir/GenomeIndex
cat $GenomeM > Genome.fa
samtools faidx Genome.fa

###Protein location
miniprot -t $Thread --gff Genome.fa $OUTdir/ProteinDB/ProteinDB.RM.fa > $OUTdir/Miniprot.gff
cd $OUTdir
cat Miniprot.gff | grep "##PAF" | awk '{print $7,$9-3000,$10+3000,"ID"NR,$2}' | sed 's/ /\t/g' > Run.list.txt

List=$OUTdir/Run.list.txt
GenomeIndex=$OUTdir/GenomeIndex/Genome.fa
ProteinDB=$OUTdir/ProteinDB/ProteinDB.RM.fa

mkdir $OUTdir/RunOUT
for Sample in $(cat $List | sed 's/\t/:/g')
do
SRA=(${Sample//\:/ })
mkdir $OUTdir/RunOUT/${SRA[3]}/
cd $OUTdir/RunOUT/${SRA[3]}/
echo "${SRA[0]}	${SRA[1]}	${SRA[2]}" > LOC.bed
bedtools getfasta -s -fi $GenomeIndex -bed LOC.bed -fo LOC.bed.fa
echo "${SRA[4]}" > Protein.ID.txt
seqtk subseq $ProteinDB Protein.ID.txt > Protein.fa
ID="G${SRA[3]}"
$Exonerate -q Protein.fa -t LOC.bed.fa --model protein2genome --showvulgar no --softmaskquery no --softmasktarget yes --showalignment no --showtargetgff yes --showcigar no --score 250 --bestn 1 --percent 50 --geneseed 250 --verbose 0 --gff3 yes > exonerate.gff
grep -v "#" exonerate.gff | sed 's/match/gene/g' | sed 's/gene_part/exon/g' | awk -F "\t" 'ARGIND==1{C=$1;S=$2;E=$2}ARGIND==2{print C"\texonerate\t"$3"\t"S+$4"\t"S+$5"\t"$6"\t"$7"\t"$8"\t"$9}' LOC.bed - | sed 's/gene00001/'$ID'/g' > exonerate.G.gff
cat exonerate.G.gff >> $OUTdir/Exonerate.gff
gth -genomic LOC.bed.fa -protein Protein.fa -intermediate -gff3out > Pudorinus.gff
grep -v "#" Pudorinus.gff | awk 'ARGIND==1{C=$1;S=$2;E=$2}ARGIND==2{print $1"\t"$2"\t"$3"\t"S+$4"\t"S+$5"\t"$6"\t"$7"\t"$8"\t"$9}' LOC.bed - | sed 's/gene1/'$ID'/g' > Pudorinus.G.gff
cat Pudorinus.G.gff >> $OUTdir/Pxonerate.gff
done



 
