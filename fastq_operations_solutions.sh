#### Download the datasets from https://zenodo.org/records/14258052

### Structure of fastq files
# fastq files are the sequencing output files files obtained after demultimplexing data from illumina seqeuncers. Most sequencing companies provide this file format for illumina data.
# It consists of four lines


@SRR15369215.126490887 #sequence identifier 
GGACCTTCTGTCATTTCACTCCTTCTGAAGTAAGGAGTGAAGTAAACACGAAGTAAACACGACAGGTTAGTCCTATTCCTTCAAGCAGGAGTACAGAAAAGAATGCAAATTCTGGGTTCTAGCCCAGCTTTTACTCCTATGGTTCTATTT #sequence
+ #seperator
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFAFJJJJJJJJJJJJJJJJFJJJJJJJJJFJJJJJJJJJFFFJJJJJJJJJJJ #base call quality scores (ASCII characters)


### Count the number of sequencing reads in the fastq files

### Solution
## 1
grep '^@' LGS1_sub_1.fq | wc -l 
# grep '^@' Finds all lines starting with @ in the file.
# ^@: Matches lines where @ appears at the beginning.
# wc -l: Counts the number of matching lines.

## 2
echo $(( $(wc -l < LGS1_sub_1.fq) / 4 ))
# echo: displays result of the arithmatic operation
# wc -1: counts the number of lines in LGS1_sub_1.fq
# dividing by 4 since one sequencing read has 4 lines of fastq file

## 3
wc -l LGS1_sub_1.fq | awk '{print $1 / 4}'

## 4
