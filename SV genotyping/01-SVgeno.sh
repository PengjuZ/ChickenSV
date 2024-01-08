#!/bin/bash
####################################################################################
###Files
#Group="10"

SampleList=./03-Cgenome/03-Genome_analysis/00-DataList/List.$Group
Index=./03-Cgenome/03-Genome_analysis/05-SVmerge/05-Graph-Con/Index
WorkDir=./03-Cgenome/03-Genome_analysis/01-SVcalling
DataDir=./03-Cgenome/00-Data/Population/Clean_data
Beagle=./00-Software/beagle.22Jul22.46e.jar
:<<Masknote
Masknote
###SoftWares
module load anaconda3/4.12.0
source activate PanCattle
fastp=./00-Software/fastp
Thread="6" #####change
####################################################################################
#####Genotyping
for Sample in $(cat $SampleList | sed 's/\t/:/g')
do
SRA=(${Sample//\:/ })
ID=${SRA[0]}
mkdir $WorkDir/$ID/SVrun
RunOUT=$WorkDir/$ID/SVrun
cat $DataDir/$ID\_clean_1.fq.gz > $RunOUT/$ID\_1.fastq.gz
cat $DataDir/$ID\_clean_2.fq.gz > $RunOUT/$ID\_2.fastq.gz
$fastp -i $RunOUT/$ID\_1.fastq.gz -I $RunOUT/$ID\_2.fastq.gz -o $RunOUT/$ID\_1.QC.fastq.gz -O $RunOUT/$ID\_2.QC.fastq.gz -j $RunOUT/$ID.fastp.json -h $RunOUT/$ID.fastp.html -q 20 -u 30 -l 75 -w $Thread
Fastq1=$RunOUT/$ID\_1.QC.fastq.gz
Fastq2=$RunOUT/$ID\_2.QC.fastq.gz
cp $Index/WChang.gg $RunOUT/$ID.gg
vg giraffe -x $Index/WChang.xg -g $RunOUT/$ID.gg -H $Index/WChang.gbwt -m $Index/WChang.min -d $Index/WChang.dist -f $Fastq1 -f $Fastq2 -t $Thread -b fast -N $ID -p > $RunOUT/$ID.gam
vg pack -Q 5 -x $Index/WChang.xg -g $RunOUT/$ID.gam -o $RunOUT/$ID.pack -t $Thread
vg call $Index/WChang.xg -r $Index/WChang.snarls -k $RunOUT/$ID.pack -s $ID -v $Index/WChang.vcf -t $Thread --bias-mode --het-bias 2,4 > $WorkDir/$ID/SV.genotype.vcf
rm -rf $WorkDir/$ID/SVrun
done
