#! /bin/bash
Evaluation=./Evaluation/
Detgaps=$Evaluation/asset-1.0.3/bin/detgaps
Ast_pb=$Evaluation/asset-1.0.3/bin/ast_pb
Ast_hic=$Evaluation/asset-1.0.3/bin/ast_hic
Acc=$Evaluation/asset-1.0.3/bin/acc
find_telomere=$Evaluation/find_telomere
threads=5

ID="Hap2"
Genome=./01-Assembly/07-ConfirmChr/$ID.fa
HIFI=./01-Assembly/02-Hybrid-trioSplit/Split/$ID.HIFI_classified.yes.fq
ONT=./01-Assembly/02-Hybrid-trioSplit/Split/$ID.ONT_classified.yes.fq
WGSAr1=./00-Data/NGS_1.QC.fq.gz
WGSAr2=./00-Data/NGS_2.QC.fq.gz
WGSr1=./01-Assembly/02-Hybrid-trioSplit/Split/$ID\_classified_524_1_.fq
WGSr2=./01-Assembly/02-Hybrid-trioSplit/Split/$ID\_classified_524_2_.fq
HICr1=./01-Assembly/02-Hybrid-trioSplit/Split/$ID.HIC_classified_590_1_.fq
HICr2=./01-Assembly/02-Hybrid-trioSplit/Split/$ID.HIC_classified_590_2_.fq
Dir=./01-Assembly/08-Assembly-evaluation/$ID
Meryl=./verkko-v1.1/lib/verkko/bin/meryl
generateDotPlot=./Evaluation/generateDotPlot

######Prepare
source /share/apps/anaconda3/bin/activate Ntools
mkdir $Dir/00-Genome
cd $Dir/00-Genome
cat $Genome > $ID.fa
Genome=$Dir/00-Genome/$ID.fa

:<<Masknote
Masknote
######Basic assembly continuity
mkdir $Dir/01-Basic
cd $Dir/01-Basic
assembly-scan $Genome > Basic.stats.txt

######Reliable blocks
mkdir $Dir/02-ReliableR
#GAP
mkdir $Dir/02-ReliableR/GAP
cd $Dir/02-ReliableR/GAP
$Detgaps $Genome > gaps.bed
samtools faidx $Genome
awk '{print $1"\t0\t"$2}' $Genome.fai > asm.bed
cat asm.bed | awk '($3-$2) > 2000 {print $1"\t0\t1000\n"$1"\t"($3-1000)"\t"$3}' > asm.ends.bed
#HiFi
mkdir $Dir/02-ReliableR/Platform
cd $Dir/02-ReliableR/Platform
minimap2 -t $threads -x map-pb -d $Genome.idx $Genome
minimap2 -x map-pb -t $threads $Genome.idx $HIFI > Genome.HiFi.paf
$Ast_pb -M 300 Genome.HiFi.paf > Genome.HiFi.bed
bedtools subtract -a $Dir/02-ReliableR/GAP/asm.bed -b Genome.HiFi.bed | bedtools merge -d 100 -i - > Genome.HiFi.low_high.bed
bedtools subtract -a Genome.HiFi.low_high.bed -b $Dir/02-ReliableR/GAP/asm.ends.bed -A > Genome.HiFi.low_high.trim1k.bed 
bedtools subtract -a $Dir/02-ReliableR/GAP/asm.bed -b Genome.HiFi.low_high.bed > Genome.HiFi.support.bed
#ONT
minimap2 -x map-pb -t $threads $Genome.idx $ONT > Genome.ONT.paf
$Ast_pb -M 300 Genome.ONT.paf > Genome.ONT.bed
bedtools subtract -a $Dir/02-ReliableR/GAP/asm.bed -b Genome.ONT.bed | bedtools merge -d 100 -i - > Genome.ONT.low_high.bed
bedtools subtract -a Genome.ONT.low_high.bed -b $Dir/02-ReliableR/GAP/asm.ends.bed -A > Genome.ONT.low_high.trim1k.bed 
bedtools subtract -a $Dir/02-ReliableR/GAP/asm.bed -b Genome.ONT.low_high.bed > Genome.ONT.support.bed
#Hi-C
bwa index $Genome
bwa mem -t $threads $Genome $HICr1 $HICr2 | samtools view -b - > Genome.HIC.bam
$Ast_hic $Dir/02-ReliableR/GAP/gaps.bed Genome.HIC.bam > Genome.HIC.bed
bedtools subtract -a $Dir/02-ReliableR/GAP/asm.bed -b Genome.HIC.bed | bedtools merge -d 100 -i - > Genome.HIC.low_high.bed
bedtools subtract -a Genome.HIC.low_high.bed -b $Dir/02-ReliableR/GAP/asm.ends.bed -A > Genome.HIC.low_high.trim1k.bed 
bedtools subtract -a $Dir/02-ReliableR/GAP/asm.bed -b Genome.HIC.low_high.bed > Genome.HIC.support.bed
#WGS
bwa mem -t $threads $Genome $WGSr1 $WGSr2 | samtools view -b - > Genome.WGS.bam
$Ast_hic $Dir/02-ReliableR/GAP/gaps.bed Genome.WGS.bam > Genome.WGS.bed
bedtools subtract -a $Dir/02-ReliableR/GAP/asm.bed -b Genome.WGS.bed | bedtools merge -d 100 -i - > Genome.WGS.low_high.bed
bedtools subtract -a Genome.WGS.low_high.bed -b $Dir/02-ReliableR/GAP/asm.ends.bed -A > Genome.WGS.low_high.trim1k.bed 
bedtools subtract -a $Dir/02-ReliableR/GAP/asm.bed -b Genome.WGS.low_high.bed > Genome.WGS.support.bed
#Reliable blocks
$Acc $Dir/02-ReliableR/GAP/gaps.bed $Dir/02-ReliableR/Platform/*.support.bed > $Dir/02-ReliableR/acc.bed

#awk '$4>1' acc.bed | bedtools merge -i - > acc.gt2.mrg.bed
# Get low support regions by merging blocks <100bp apart
#bedtools subtract -a asm.bed -b acc.gt2.mrg.bed | bedtools merge -d 100 -i - > low_support.bed
# Get the final support region as reliable blocks
#bedtools subtract -a asm.bed -b low_support.bed > reliable.bed
# Exclude low supports in <1kb scaffold boundaries for excluding end-scaffold effects
#bedtools subtract -A -a low_support.bed -b asm.ends.bed > low_support.trim1k.bed

######Telomere
mkdir $Dir/03-Telomere
cd $Dir/03-Telomere
$find_telomere $Genome > Genome.telomere
sdust $Genome > Genome.sdust
source /share/apps/anaconda3/bin/activate nextDenovo
java -cp $Evaluation/telomere8.jar SizeFasta $Genome > Genome.lens
java -cp $Evaluation/telomere8.jar FindTelomereBreaks Genome.lens Genome.sdust Genome.telomere > Genome.breaks
# Lowering threshold to 0.10 (10%) from the initial 0.40 (40%)
threshold=0.4
java -cp $Evaluation/telomere8.jar FindTelomereWindows Genome.telomere 99.9 $threshold > Genome.windows.$threshold
#Merge telomere motifs in 100bp
cat Genome.windows.$threshold | awk '{print $2"\t"$(NF-2)"\t"$(NF-1)}' | sed 's/>//g' | bedtools merge -d 100  > Genome.windows.$threshold.bed
#Find those at end of scaffolds, within < $ends
ends=10000
cat Genome.lens | awk -v ends=$ends '{if ($2>(ends*2)) {print $1"\t0\t"ends"\n"$1"\t"($NF-ends)"\t"$NF} else {print $1"\t0\t"$NF}}' > Genome.ends.bed
bedtools intersect -wa -a Genome.windows.$threshold.bed -b Genome.ends.bed > Genome.windows.$threshold.$ends.ends.bed

######Merqury
mkdir $Dir/04-Merqury
cd $Dir/04-Merqury
source /share/apps/anaconda3/bin/activate merqury
$Meryl k=21 count output WGS.R1.meryl $WGSAr1
$Meryl k=21 count output WGS.R2.meryl $WGSAr2
meryl union-sum output Genome.meryl WGS.*.meryl
merqury.sh Genome.meryl $Genome OUT

######BUSCO
mkdir $Dir/05-Busco
cd $Dir/05-Busco
source /share/apps/anaconda3/bin/activate Busco
#busco -m genome -i $Genome -o OUTPUT -l ./Busco/eukaryota_odb10 --offline -f -c 32
busco -m genome -i $Genome -o OUTPUT -l ./Busco/aves_odb10 --offline -f -c $threads

######genome-to-genome alignment
mkdir $Dir/06-G2G
cd $Dir/06-G2G
ref=./00-Chicken-NCBI/GCF_016699485.2/Genome.fa
###mashmap
# --pi 90 ; -s 2000
mashmap -r $ref -q $Genome -t $threads -o G2G.map --filter_mode one-to-one 
perl $generateDotPlot png large G2G.map
cat G2G.map | awk -F " " -v name=$name '{print $6"\t"$8"\t"$9+1"\t"$1"\t"$3"\t"$4+1"\t"$1":"$3"-"$4+1":"name"\t"$NF"\t"$5"\t"$7"\t"$2}' > G2G.map.bed
###nucmer
nucmer -p nucmer --maxmatch -l 100 -c 500 $ref $Genome
dnadiff -p nucmer.dnadiff -d nucmer.delta
mummerplot -p nucmer.dnadiff.plot --layout --fat -t png nucmer.dnadiff.1delta






