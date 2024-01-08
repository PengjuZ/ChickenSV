#!/bin/sh
#wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.9/Rfam.cm.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.9/Rfam.clanin
#conda create -n NonCoding
#conda install -c bioconda infernal
#conda install -c bioconda trnascan-se

GenomeUM=$1
Workdir=$2
Thread=$3
RfamDB=./03-Cgenome/00-Chicken-NCBI/Rfam/
source /share/apps/anaconda3/bin/activate NonCoding

###infernal
mkdir $Workdir/Infernal
cd $Workdir/Infernal
esl-seqstat $GenomeUM | grep "residues" | awk '{print $4}' > Seqstat.txt
less $RfamDB/Rfam.cm | grep 'NAME'| sort | wc -l > Rfamstat.txt
paste Seqstat.txt Rfamstat.txt | awk '{print int($1*$2*2/1000000)}' - > Input.txt
var=$(cat Input.txt)
echo ${var}
cmscan -Z ${var} --cut_ga --rfam --nohmmonly --tblout OUT.tblout --fmt 2 --cpu $Thread --clanin $RfamDB/Rfam.clanin $RfamDB/Rfam.cm $GenomeUM > OUT.cmscan
awk 'BEGIN{OFS="\t";}{if(FNR==1) print "target_name\taccession\tquery_name\tquery_start\tquery_end\tstrand\tscore\tEvalue"; if(FNR>2 && $20!="=" && $0!~/^#/) print $2,$3,$4,$10,$11,$12,$17,$18; }' OUT.tblout > OUT.tblout.final.xls
###tRNAscan-SE
mkdir $Workdir/tRNAscan
cd $Workdir/tRNAscan
tRNAscan-SE -o tRNA.out -f tRNA.ss -m tRNA.stats $GenomeUM
