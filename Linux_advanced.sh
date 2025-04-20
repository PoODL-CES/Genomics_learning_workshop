#### Download the datasets from https://zenodo.org/records/14258052
#fastq reads are in .gz format
# Original genome files are in .fna file (fna is fasta)
#fna.amb: stores information about the presence of ambiguous bases (e.g., N, or other non-standard nucleotides) in the reference genome.
#fna.ann:  stores information about the reference genome, including sequence names, lengths, and other metadata.
#fna.bwt: contains the Burrows-Wheeler transformed sequence, a compressed representation of the reference genome.
#fna.pac: contains packed sequence data
#fna.sa: suffix array index, used to locate positions of sequences within the reference genome.

### Structure of fastq files
# fastq files are the sequencing output files obtained after demultimplexing data from illumina seqeuncers. Most sequencing companies provide this file format for illumina data in ".gz" extension.It is primarily used for reducing the size of a single file, often on UNIX-like systems.
# It consists of four lines


@SRR15369215.126490887 #sequence identifier 
GGACCTTCTGTCATTTCACTCCTTCTGAAGTAAGGAGTGAAGTAAACACGAAGTAAACACGACAGGTTAGTCCTATTCCTTCAAGCAGGAGTACAGAAAAGAATGCAAATTCTGGGTTCTAGCCCAGCTTTTACTCCTATGGTTCTATTT #sequence
+ #seperator
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFAFJJJJJJJJJJJJJJJJFJJJJJJJJJFJJJJJJJJJFFFJJJJJJJJJJJ #base call quality scores (ASCII characters)


#### Count the number of sequencing reads in the fastq files

#### How many reads are shorter than 150bp?

#### How many reads are longer than 150bp?

#### How many reads are not equal to 150bp?

#### How many reads are contaminated with illumination adapters
#### The adapter sequence is CTGTCTCTTATACACATCT
