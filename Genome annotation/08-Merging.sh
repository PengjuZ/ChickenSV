#!/bin/sh
#conda create -n Merging
#conda install -c bioconda evidencemodeler
#conda install -c bioconda bedtools

GenomeUM=$1
Workdir=$2
Rundir=$3
Thread=$4
Refanno=$5
RNAdb=$6
ID=$7

ToolsDir=./.conda/envs/Merging/opt/evidencemodeler-1.1.1/EvmUtils/
Gffrename=./00-Software/gffrename.py

Augustus=./03-Cgenome/02-Annotation/$ID/04-Abinitio/Augustus/Augustus.gff
Genemark=./03-Cgenome/02-Annotation/$ID/04-Abinitio/GeneMark/Genemark.gff
Braker2=./03-Cgenome/02-Annotation/$ID/05-RNA-seq/Alignments/OUT/augustus.hints_utr.gff3
Transdecoder=./03-Cgenome/02-Annotation/$ID/05-RNA-seq/Assembly/OUT/transcripts.fasta.transdecoder.genome.gff3
Exonerate=./03-Cgenome/02-Annotation/$ID/06-Homology/Exonerate.gff
Pxonerate=./03-Cgenome/02-Annotation/$ID/06-Homology/Pxonerate.gff

##############Weights

source /share/apps/anaconda3/bin/activate Merging
mkdir $Rundir/Partition
cd $Rundir/Partition
cat $Augustus $Genemark > Gene_predictions.gff3
cat $Braker2 | sed 's/AUGUSTUS/OTHER_PREDICTION/g' >> Gene_predictions.gff3
cat $Transdecoder | sed 's/transdecoder/OTHER_PREDICTION/g' >> Gene_predictions.gff3
cat $Exonerate $Pxonerate | awk '{if($3 == "gene"){print}else{T[$3][$9]+=1;R="ID="$3T[$3][$9]";"$9;$9=R;print}}' | sed 's/ /\t/g' > Proteins.gff3


$ToolsDir/partition_EVM_inputs.pl \
--genome $GenomeUM \
--gene_predictions Gene_predictions.gff3 \
--protein_alignments Proteins.gff3 \
--segmentSize 1000000 \
--overlapSize 100000 \
--partition_listing partitions_list.out


$ToolsDir/write_EVM_commands.pl \
--genome $GenomeUM \
--weights $Rundir/weights.txt \
--min_intron_length 10 \
--gene_predictions Gene_predictions.gff3 \
--protein_alignments Proteins.gff3 \
--output_file_name evm.out \
--partitions partitions_list.out > Commands.list


parallel --jobs $Thread < $Rundir/Partition/Commands.list


$ToolsDir/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out


$ToolsDir/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm.out --genome $GenomeUM


find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff


