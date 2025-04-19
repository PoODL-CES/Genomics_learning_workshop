####Task 1: Map trimmed reads using bwa mem
## The datasets from https://zenodo.org/records/14258052 contain the reference genome file and its index files

### Install bwa
conda create -n bwa -c bioconda bwa

## Activate the conda environment
conda activate bwa

#aligning reads with bwa mem (burrows Burrows-Wheeler Aligner)
bwa mem GCA_021130815.1_PanTigT.MC.v3_genomic.fna BEN_NW10_sub_1_val_1.fq.gz BEN_NW10_sub_2_val_2.fq.gz > BEN_NW_10_aligned_reads.sam
#bwa mem: runs the "mem" algorithm of BWA. It is optimum for 70bp-1Mbp reads, and commonly used for Illumina short-read data.
#GCA_021130815.1_PanTigT.MC.v3_genomic.fna: reference genome in FASTA format. The bwa index file should also be present.
#../fq_files/BEN_NW10_sub_1_val_1.fq.gz ../fq_files/BEN_NW10_sub_2_val_2.fq.gz: These are the paired-end FASTQ files.
#> BEN_NW_10_aligned_reads.sam: SAM file alignment output. It contains detailed alignment information for each read.


# Mapping all reads to reference genome in single step

for file1 in *_sub_1_val_1.fq.gz; do
    file2=${file1/_sub_1_val_1.fq.gz/_sub_2_val_2.fq.gz}
    sample_name=$(basename "$file1" _sub_1_val_1.fq.gz)
    
    bwa mem GCA_021130815.1_PanTigT.MC.v3_genomic.fna "$file1" "$file2" > "${sample_name}_aligned_reads.sam"
done


## Deactivate the conda environment
conda deactivate

#Convert sam to bam 
# Install samtools
conda create -n samtools -c bioconda samtools

## EXPLAIN THE SAMTOOLS OPTIONS USED. WHY IS THIS STEP NECESSARY?
samtools view -S -b BEN_NW_10_aligned_reads.sam | samtools sort -o BEN_NW_10_sorted_reads.bam
### view: This is the subcommand in samtools used to convert, filter, or view alignment files.
## -S: Specifies that the input file is in SAM (Sequence AlignmentMap) format.
## -b: Specifies that the output file should be in BAM (Binary Alignment Map) format, which is a compressed version of the SAM format.

## We use BAM files because these files are much smaller in size compared to SAM files, saving storage space.
## They are compressed and indexed, which allows faster access and processing during downstream analysis.

## Convert all sam files to bam and sort the reads

for file in *.sam; do 
    samtools view -S -b "$file" | samtools sort -o "${file%.sam}_sorted.bam"
done

## sort: This samtools subcommand is used to sort the alignment data based on genomic coordinates.
# Sorting the reads improves the efficiency of downstream analysis like variant calling etc. by giving access to specific regions of the genome.

## Deactivate the conda environment
conda deactivate 

##### MARKDUPLICATES USING GATK4
### Install and activate gatk4
conda create -n gatk4 -c bioconda gatk4
conda activate gatk4

### Mark and remove duplicates
gatk MarkDuplicates -I BEN_NW_10_sorted_reads.bam -O BEN_NW_10_deduplicated.bam -M BEN_NW_10_duplication_metrics.txt --REMOVE_DUPLICATES true
# Markduplicates is necessary for marking duplicate reads which arise during pcr amplification step

## Marking and removing duplicates all at once
## we are using the parallel command here because of small group size and small files. 
## However this is memory intensive and computers can crash if run by big groups on heavy files
parallel 'gatk MarkDuplicates -I {} -O {.}_deduplicated.bam -M {.}_duplication_metrics.txt --REMOVE_DUPLICATES true' ::: *_sorted.bam

### INDEX AFTER MARKDUPLICATES
samtools index BEN_NW_10_deduplicated.bam
# indexing allows quick access to specific genomic regions and improve performance of downstream analysis tools

## Indexing the deduplicated files all at once
parallel 'samtools index {}' ::: *_deduplicated.bam

## for statistics file
parallel 'samtools stats {} > {.}_stats.txt' ::: *_deduplicated.bam
## This will generate a statistics file which will have information about the number of reads that mapped to the reference geome, number of unmapped reads etc.

#### for estimating sequencing statistics like coverage per chromosome/scaffold
conda install -c bioconda qualimap

qualimap bamqc -bam BEN_NW12_aligned_reads_sorted_deduplicated.bam -outdir qualimap_results -outformat HTML

#-bam : to input bam file
#-outdir : Directory for results
#-outformat HTML :  Output in HTML 

cd qualimap_results
cat genome_results.txt

# Look for a section Chromosome-wise coverage” or similar
