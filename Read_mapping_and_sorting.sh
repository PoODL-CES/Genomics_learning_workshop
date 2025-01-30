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
### view: This is the subcommand in samtools used to convert, filter, or view alignment files.
## -S: Specifies that the input file is in SAM (Sequence AlignmentMap) format.
## -b: Specifies that the output file should be in BAM (Binary Alignment Map) format, which is a compressed version of the SAM format.

## We use BAM files because these files are much smaller in size compared to SAM files, saving storage space.
## They are compressed and indexed, which allows faster access and processing during downstream analysis.

## EXPLAIN THE SAMTOOLS OPTIONS USED. WHY IS THIS STEP NECESSARY?
samtools sort BEN_NW_10_aligned_reads.bam -o BEN_NW_10_sorted_reads.bam
## sort: This samtools subcommand is used to sort the alignment data based on genomic coordinates.
# Sorting the reads improves the efficiency of downstream analysis like variant calling etc. by giving access to specific regions of the genome.

## Deactivate the conda environment

### MARKDUPLICATES USING PICARD

### INDEX AFTER MARKDUPLICATES USING PICARD

samtools index BEN_NW_10_sorted_reads.bam
samtools flagstat BEN_NW_10_sorted_reads.bam
