#### Map reads using bwa mem
### Install bwa
conda create -n bwa -c bioconda bwa

## Activate the conda environment
conda activate bwa

bwa mem GCA_021130815.1_PanTigT.MC.v3_genomic.fna ../fq_files/BEN_NW10_sub_1_val_1.fq.gz ../fq_files/BEN_NW10_sub_2_val_2.fq.gz > BEN_NW_10_aligned_reads.sam

## Deactivate the conda environment

### Install samtools
FILL UP THE INSTALLATION


## EXPLAIN THE SAMTOOLS OPTIONS USED. WHY IS THIS STEP NECESSARY?
samtools view -S -b BEN_NW_10_aligned_reads.sam > BEN_NW_10_aligned_reads.bam


## EXPLAIN THE SAMTOOLS OPTIONS USED. WHY IS THIS STEP NECESSARY?
samtools sort BEN_NW_10_aligned_reads.bam -o BEN_NW_10_sorted_reads.bam

## Deactivate the conda environment

### MARKDUPLICATES USING PICARD





samtools index BEN_NW_10_sorted_reads.bam
samtools flagstat BEN_NW_10_sorted_reads.bam
