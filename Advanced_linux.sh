#### Download the datasets from https://zenodo.org/records/14258052

### Structure of fastq files
# fastq files are the sequencing output files files obtained after demultimplexing data from illumina seqeuncers. Most sequencing companies provide this file format for illumina data.
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
