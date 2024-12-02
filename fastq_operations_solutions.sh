#### Download the datasets from https://zenodo.org/records/14258052

### Structure of fastq files
# fastq files are the sequencing output files files obtained after demultimplexing data from illumina seqeuncers. Most sequencing companies provide this file format for illumina data.
# It consists of four lines


@SRR15369215.126490887 #sequence identifier 
GGACCTTCTGTCATTTCACTCCTTCTGAAGTAAGGAGTGAAGTAAACACGAAGTAAACACGACAGGTTAGTCCTATTCCTTCAAGCAGGAGTACAGAAAAGAATGCAAATTCTGGGTTCTAGCCCAGCTTTTACTCCTATGGTTCTATTT #sequence
+ #seperator
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFAFJJJJJJJJJJJJJJJJFJJJJJJJJJFJJJJJJJJJFFFJJJJJJJJJJJ #base call quality scores (ASCII characters)


#### Count the number of sequencing reads in the fastq files

### Solution
## 1
less LGS1_sub_1.fq.gz | grep '^@' | wc -l 
# grep '^@' Finds all lines starting with @ in the file.
# ^@: Matches lines where @ appears at the beginning.
# wc -l: Counts the number of matching lines.

## 2
echo $(( $(zcat LGS1_sub_1.fq.gz | wc -l) / 4 ))
#zcat: to decompress .gz file
#wc-1: counts the number of lines in decompressed output
#divsion by 4 to calculate total number of reads (since each fastq entry will have 4 lines)
#echo to show the result of the calculation

## 3
less LGS1_sub_1.fq.gz | wc -l | awk '{print $1 / 4}'

#### How many reads are shorter than 150bp?
