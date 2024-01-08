#!/bin/bash
####################################################################################
#conda create -n GenomeE -c bioconda mashtree
#conda create -n Pggb
#conda install -c bioconda wfmash
#conda install -c bioconda seqwish
#conda install -c bioconda smoothxg
#conda install -c bioconda odgi
#conda install -c bioconda gfaffix
#conda install -c bioconda bcftools
#conda install -c bioconda vg
#conda install -c bioconda multiqc
#conda install -c conda-forge pigz
#conda install -c bioconda vcflib
#conda install -c bioconda pggb



Dir=./08-ChiGtex
Workdir=./08-ChiGtex/03-PanCon


source /share/apps/anaconda3/bin/activate Ntools
:<<Masknote
Masknote
##################################
###01-Merge
mkdir $Workdir/01-Merge/
cd $Workdir/01-Merge/
ID="G01 G02 G03 G04 G05 G06 G07 G08 G09 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20 G21 G22 G23 G24 G25 G26 G27 G28 G29 G30 G31"

for i in $ID
do
awk '{if(substr($1,0,1) == ">"){print ">"I"#1#"substr($1,2)}else{print}}' I=$i $Dir/00-Data/$i.fa | seqkit seq -u - >> Merge.fa
done
bgzip -@ 64 Merge.fa
samtools faidx Merge.fa.gz

##################################
###02-Split

ID="G01 G02 G03 G04 G05 G06 G07 G08 G09 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20 G21 G22 G23 G24 G25 G26 G27 G28 G29 G30 G31"
mkdir $Workdir/02-Split/
cd $Workdir/02-Split/
for i in $ID
do
echo $i
awk '{if($10 >= 90){S[$1][$6]+=$4-$3;print $1"\t"$6"\t"S[$1][$6]"\t"$2}}' $Dir/02-Evaluation/$i.G2G.map | sort -k1,1 -k2,2 -k3,3nr | awk '{if(M[$1][$2]!=1){print $1,$2,$3,$4,$3/$4;M[$1][$2]=1}}' > $i.link.txt
sed 's/ /\t/g' $i.link.txt | sort -k1,1 -k3nr,3 | awk '{if(M[$1]!=1){print I"#1#"$1,"G01#1#"$2,$3,$4,$5 ;M[$1]=1}}' I=$i | sed 's/ /\t/g' > $i.link.yes.txt
done

##################################
###03-Partitioning
source /share/apps/anaconda3/bin/activate Pggb
mkdir $Workdir/03-Partitioning/
cd $Workdir/03-Partitioning/
wfmash $Workdir/01-Merge/Merge.fa.gz -p 98 -n 31 -t 50 -s 10000 -m > Merge.mapping.paf
cat $Workdir/02-Split/*.link.yes.txt | awk 'ARGIND==1{M[$1][$2]=1;M[$2][$1]=1;S[$1]=1}ARGIND==2{if(M[$1][$6]==1){print}}' - Merge.mapping.paf > Merge.filter.paf
paf2net.py -p Merge.filter.paf
net2communities.py \
    -e Merge.filter.paf.edges.list.txt \
    -w Merge.filter.paf.edges.weights.txt \
    -n Merge.filter.paf.vertices.id2name.txt \
    --plot
seq 0 40 | while read i; do
    samtools faidx $Workdir/01-Merge/Merge.fa.gz $(cat Merge.filter.paf.edges.weights.txt.community.$i.txt) | \
    bgzip -@ 64 -c > Merge.filter.community.$i.fa.gz
    samtools faidx Merge.filter.community.$i.fa.gz
done

##################################
###04-GraphBuilding
mkdir $Workdir/04-GraphBuilding/
seq 0 40 | while read i; do
mkdir $Workdir/04-GraphBuilding/$i/
cd $Workdir/04-GraphBuilding/$i/
Number=$(sed 's/#/\t/g' $Workdir/03-Partitioning/Merge.filter.community.$i.fa.gz.fai | cut -f 1 | uniq | wc -l)
echo $Number
pggb -i $Workdir/03-Partitioning/Merge.filter.community.$i.fa.gz -o output -n $Number -t 64 -p 95 -s 10000 -T 20 --poa-params 1,9,16,2,41,1
cp $Workdir/04-GraphBuilding/$i/output/*.gfa Pan.gfa
vg deconstruct -P 'G04' -H '#' -e -a -t 64 Pan.gfa > Pan.G04.vcf
done


