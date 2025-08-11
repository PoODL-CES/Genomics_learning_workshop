####Task 1: Map trimmed reads using bwa mem (Trimming-cleaning the raw sequences; Mapping-aligning to the reference genome)
## The datasets from https://zenodo.org/records/14258052 contain the reference genome file and its index files

### Install bwa
conda create -n bwa -c bioconda bwa
BWA (burrows Burrows-Wheeler)- Aligner tool used to align DNA sequencing reads to a reference genome.

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
#for file1 in *_sub_1_val_1.fq.gz; do: looks through all the read 1(forward) FASTQ files.
#file2=${file1/_sub_1_val_1.fq.gz/_sub_2_val_2.fq.gz}: constructs the corresponding read 2 (reverse) file name by replacing _sub_1_val_1.fq.gz with _sub_2_val_2.fq.gz in file1.
#sample_name=$(basename "$file1" _sub_1_val_1.fq.gz): removes the _sub_1_val_1.fq.gz part from the filename, leaving just the sample ID (e.g., BEN_NW10).
# bwa mem GCA_021130815.1_PanTigT.MC.v3_genomic.fna "$file1" "$file2" > "${sample_name}_aligned_reads.sam": runs the actual alignment

## Deactivate the conda environment
conda deactivate

#Task 2: Convert sam to bam 
# Install samtools
conda create -n samtools -c bioconda samtools
conda activate samtools
## EXPLAIN THE SAMTOOLS OPTIONS USED. WHY IS THIS STEP NECESSARY?
samtools view -S -b BEN_NW_10_aligned_reads.sam | samtools sort -o BEN_NW_10_sorted_reads.bam
### view: This is the subcommand in samtools used to convert, filter, or view alignment files.
## -S: Specifies that the input file is in SAM (Sequence AlignmentMap) format.
## -b: Specifies that the output file should be in BAM (Binary Alignment Map) format, which is a compressed version of the SAM format.

## We use BAM files because these files are much smaller in size compared to SAM files, saving storage space.
## They are compressed and indexed, which allows faster access and processing during downstream analysis.

## (Task 3) Convert all sam files to bam and sort the reads

for file in *.sam; do 
    samtools view -S -b "$file" | samtools sort -o "${file%.sam}_sorted.bam"
done

#samtools view: samtools view converts SAM file to BAM format.
#-S -b "$file": -s tells samtools that the input is in SAM format. -b tells that the output should be in BAM format. "$file" is the name of the current .sam file in the loop.
#samtools sort: sorts the BAM file by genomic coordinates.
#"${file%.sam}: strips off the .sam extension
#then _sorted.bam is added.
## sort: This samtools subcommand is used to sort the alignment data based on genomic coordinates.
#Sorting necessary for BAM file indexing, downstream analysis, visualization, duplicate marking.

## Deactivate the conda environment
conda deactivate 

## Task 4: MARKDUPLICATES USING GATK4
### Install and activate gatk4
conda create -n gatk4 -c bioconda gatk4
conda activate gatk4

### Mark and remove duplicates
gatk MarkDuplicates -I BEN_NW_10_sorted_reads.bam -O BEN_NW_10_deduplicated.bam -M BEN_NW_10_duplication_metrics.txt --REMOVE_DUPLICATES true
# Markduplicates is necessary for marking duplicate reads which arise during pcr amplification step

## Marking and removing duplicates all at once
## we are using the parallel command here because of small group size and small files. 
## However this is memory intensive and computers can crash if run by big groups on heavy files
conda create -n parallel -c bioconda parallel
conda activate parallel
parallel 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate bwa && gatk MarkDuplicates -I {} -O {.}_deduplicated.bam -M {.}_duplication_metrics.txt --REMOVE_DUPLICATES true' ::: *_sorted.bam
#parallel '...' ::: *_sorted.bam: takes each sorted bam file, replaces {} with the file name, replaces {.} with the filename without the .bam extension, run the specified command in parallel for each file
#'source ~/miniconda3/etc/profile.d/conda.sh':Loads Conda’s environment setup script. Without this, conda activate won’t work in a non-interactive shell like the ones parallel launches.
#gatk MarkDuplicates: 
#-I {} — input is the BAM file.
#-O {.}_deduplicated.bam — output file with duplicates removed.
#-M {.}_duplication_metrics.txt — metrics file.
#--REMOVE_DUPLICATES true — tells GATK to physically remove duplicates from the output BAM file, not just mark them.

#or

for file in *_sorted.bam; do
    base=${file%_sorted.bam}
    gatk MarkDuplicates \
        -I "$file" \
        -O "${base}_deduplicated.bam" \
        -M "${base}_duplication_metrics.txt" \
        --REMOVE_DUPLICATES true
done

#base=${file%_sorted.bam}: ${file%_sorted.bam} removes _sorted.bam from the filename.
#-I "$file": takes the input *_sorted.bam file
#-O "${base}_deduplicated.bam" \: outputs the deduplicated BAM (with duplicates removed).
#-M "${base}_duplication_metrics.txt": writes duplication statistics to this file.
#--REMOVE_DUPLICATES true: removes the duplicate reads instead of just marking them.

#files would be viewed as BEN_NW_10_sorted_reads.bam 

##Task 5: INDEX AFTER MARKDUPLICATES
samtools index BEN_NW_10_deduplicated.bam
# indexing allows quick access to specific genomic regions and improve performance of downstream analysis tools

## Indexing the deduplicated files all at once

for file in *_deduplicated.bam; do
samtools index file
done


#Task 6: for estimating sequencing statistics like coverage per chromosome/scaffold
conda create -n qualimap -c bioconda qualimap
conda activate qualimap

qualimap bamqc -bam BEN_NW12_deduplicated.bam -outdir qualimap_results -outformat HTML

#for bulk statistics

qualimap bamqc -bam *_deduplicated.bam -outdir qualimap_results -outformat HTML

#-bam : to input bam file
#-outdir : Directory for results
#-outformat HTML :  Output in HTML 

cd qualimap_results
cat genome_results.txt

# Look for a section Chromosome-wise coverage” or similar
