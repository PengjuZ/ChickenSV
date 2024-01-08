#! /bin/bash
##############################################
### Authors: Ivan Sovic and Ann Mc Cartney ###
### Date: September 2021                   ###
##############################################
### Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies
### Dependencies: meryl, minimap2, merfin, bcftools, racon, winnowmap, falconc (pbipa)
#conda create -n Polish
#source activate Polish
#conda install -c bioconda meryl
#conda install -c bioconda minimap2
#conda install -c bioconda merfin
#conda install -c bioconda bcftools
#conda install -c bioconda samtools
#conda install -c bioconda winnowmap
#conda install -c bioconda pb-falconc
:<<Masknote
mkdir ./01-Assembly/06-Polishing-Consensus/Hap2
mkdir ./01-Assembly/06-Polishing-Consensus/Hap2/HIFI
cd ./01-Assembly/06-Polishing-Consensus/Hap2/HIFI

source /share/apps/anaconda3/bin/activate Polish
# Input parameters.
CMD_OPT_NUM_THREADS=8
CMD_OPT_ITERATIONS=3
CMD_OPT_IN_DRAFT_FASTA=./01-Assembly/05-Gapfilling-trio/Hap2/07-Purging/Cgenome.fa
CMD_OPT_IN_READS=./01-Assembly/02-Hybrid-trioSplit/Split/Hap2.HIFI_classified.yes.fq
CMD_OPT_IN_READMERS=./01-Assembly/02-Hybrid-trioHIC/meryl/hap2.meryl
CMD_OPT_OUT_PREFIX=./01-Assembly/06-Polishing-Consensus/Hap2/HIFI/Hap

# Dependencies.
RACON=/share/home/zju_zhaopj/00-Software/racon/racon/build/bin/racon
WINNOWMAP=winnowmap
FALCONC=falconc
MERYL=meryl
MERFIN=merfin
BCFTOOLS=bcftools

function run_one_iteration {
    local out_prefix=$1; shift
    local in_draft=$1; shift
    local threads=$1; shift
    local in_dataset=$1; shift
    local in_readmers=$1; shift         # IlluminaPCRfree.k21.gt1.meryl

    # Get the absolute paths.
    mkdir -p $(dirname ${out_prefix})
    out_prefix=$(realpath ${out_prefix})

    # Generate repeitive 15-mers to downweight.
    local out_winnowmap_bam=${out_prefix}.winnowmap.bam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_winnowmap_bam}.count15k.memtime \
    ${MERYL} count k=15 ${in_draft} output merylDB
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_winnowmap_bam}.print15k.memtime \
    ${MERYL} print greater-than distinct=0.9998 merylDB > ${out_winnowmap_bam}.repetitive_k15.txt

    # Map the reads using Winnowmap.
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_winnowmap_bam}.memtime \
    ${WINNOWMAP} --MD -W ${out_winnowmap_bam}.repetitive_k15.txt -ax map-pb -t ${threads} ${in_draft} ${in_dataset} | samtools view -Sb > ${out_winnowmap_bam}

    # Sort the BAM file.
    local out_winnowmap_sorted_bam=${out_prefix}.winnowmap.sorted.bam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_winnowmap_sorted_bam}.memtime \
    samtools sort --threads ${threads} -o ${out_winnowmap_sorted_bam} ${out_winnowmap_bam}

    # Filtering the BAM file.
    local out_falconc_sam=${out_prefix}.falconc.sam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_falconc_sam}.memtime \
    falconc bam-filter-clipped -t -F 0x104 --input-fn ${out_winnowmap_sorted_bam} --output-fn ${out_falconc_sam} --output-count-fn ${out_falconc_sam}.filtered_aln_count.txt 2>&1 | tee ${out_falconc_sam}.falconc.log

    # Polish using Racon.
    local out_racon_fasta=${out_prefix}.racon.fasta
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_racon_fasta}.memtime \
    ${RACON} -t ${threads} ${in_dataset} ${out_falconc_sam} ${in_draft} -L ${out_racon_fasta} -S > ${out_racon_fasta}

    # Generate the Meryl database.
    local out_meryl=${out_prefix}.racon.meryl
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_meryl}.memtime \
    ${MERYL} count k=21 ${in_draft} output ${out_meryl}

    # Run Merfin.
    local out_merfin=${out_prefix}.racon.merfin
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_merfin}.memtime \
    ${MERFIN} -polish -sequence ${in_draft} -seqmers ${out_meryl} -readmers ${in_readmers} -peak 23.8 -vcf ${out_racon_fasta}.vcf -output ${out_merfin} -threads ${threads}

    # Call Consensus
    local out_consensus=${out_prefix}.consensus.fasta
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_consensus}.bcftools_view.memtime \
    ${BCFTOOLS} view -Oz ${out_merfin}.polish.vcf > ${out_merfin}.polish.vcf.gz
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_consensus}.bcftools_index.memtime \
    ${BCFTOOLS} index ${out_merfin}.polish.vcf.gz
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_consensus}.bcftools_consensus.memtime \
    ${BCFTOOLS} consensus ${out_merfin}.polish.vcf.gz -f ${in_draft} -H 1 > ${out_consensus}
}

function run_all {
    local out_prefix=$1; shift
    local threads=$1; shift
    local iterations=$1; shift
    local in_draft=$1; shift
    local in_reads=$1; shift
    local in_readmers=$1; shift         # IlluminaPCRfree.k21.gt1.meryl

    mkdir -p $(dirname ${out_prefix})
    cp ${in_draft} ${out_prefix}.iter_0.consensus.fasta
    for ((i = 0; i < ${iterations} ; i++)); do
        next_i=$(( i + 1 ))
        run_one_iteration ${out_prefix}.iter_${next_i} ${out_prefix}.iter_${i}.consensus.fasta ${threads} ${in_reads} ${in_readmers}
    done
}

run_all ${CMD_OPT_OUT_PREFIX} ${CMD_OPT_NUM_THREADS} ${CMD_OPT_ITERATIONS} ${CMD_OPT_IN_DRAFT_FASTA} ${CMD_OPT_IN_READS} ${CMD_OPT_IN_READMERS}

Masknote

##########WGS-freebayes
mkdir ./01-Assembly/06-Polishing-Consensus/Hap2/WGS
cd ./01-Assembly/06-Polishing-Consensus/Hap2/WGS
source /share/apps/anaconda3/bin/activate nextDenovo

read1=./01-Assembly/02-Hybrid-trioSplit/Split/Hap2_classified_524_1_.fq
read2=./01-Assembly/02-Hybrid-trioSplit/Split/Hap2_classified_524_2_.fq
Hap=./01-Assembly/06-Polishing-Consensus/Hap2/HIFI/Hap.iter_3.consensus.fasta
threads=20

cat $Hap > Hap.fa
input=Hap.fa
bwa index ${input}
samtools faidx ${input}
bwa mem -t ${threads} ${input} ${read1} ${read2} | samtools view --threads ${threads} -b - | samtools sort --threads ${threads} -O BAM -o aligned.bam - 
samtools index aligned.bam
freebayes --bam aligned.bam -f $input > aligned.vcf
bgzip aligned.vcf
tabix -fp vcf aligned.vcf.gz
bcftools view --no-version -Ou aligned.vcf.gz > aligned.bcf
bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=$threads aligned.bcf > aligned.changes.vcf.gz
bcftools index aligned.changes.vcf.gz
bcftools consensus -Hla -f ${input} aligned.changes.vcf.gz > Genome.fb.fasta

