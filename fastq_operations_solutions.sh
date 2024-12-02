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

### Solution
## 1
zcat BEN_CI16_sub_1.fq.gz | awk 'NR % 4 == 2 {print length($0)}' | awk '$1<150' | wc -l
## NR % 4 == 2: This will check every 2nd line of the each four line of FASTAQ and length($0): This will returns the length of the sequence.

#### How many reads are longer than 150bp?
### Solution


#### How many reads are not equal to 150bp?

### Solution
## 1
zcat file.fastq.gz | awk 'NR % 4 == 2 {if (length($0) != 150) print length($0)}' | wc -l
#if (length($0) != 150): Checks if the length of the sequence is not equal to 150.
#print length($0): Prints the length of sequences that are not 150.

# The lenght of sequence more than 150
zcat BEN_CI16_sub_1.fq.gz | awk 'NR % 4 == 2 {print length($0)}' | awk '$1>150' | wc -l 
# Similar to solution 1 of this part.

