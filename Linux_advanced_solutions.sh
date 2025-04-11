#### Download the datasets from https://zenodo.org/records/14258052

### Structure of fastq files
# fastq files are the sequencing output files files obtained after demultimplexing data from illumina seqeuncers. Most sequencing companies provide this file format for illumina data.
# It consists of four lines


@SRR15369215.126490887 #sequence identifier 
GGACCTTCTGTCATTTCACTCCTTCTGAAGTAAGGAGTGAAGTAAACACGAAGTAAACACGACAGGTTAGTCCTATTCCTTCAAGCAGGAGTACAGAAAAGAATGCAAATTCTGGGTTCTAGCCCAGCTTTTACTCCTATGGTTCTATTT #sequence
+ #seperator
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFAFJJJJJJJJJJJJJJJJFJJJJJJJJJFJJJJJJJJJFFFJJJJJJJJJJJ #base call quality scores (ASCII characters)

#After 100-150 base pairs, the accuracy of base calling decreases. Trimming is usually done to retain only high-quality bases.

# Solution

##1) Number of sequencing reads in the fastq files
less LGS1_sub_1.fq.gz | grep '^@' | wc -l 
# grep '^@' Finds all lines starting with @ in the file.
# ^@: Matches lines where @ appears at the beginning.
# wc -l: Counts the number of matching lines.

#or

#1) 
echo $(( $(zcat LGS1_sub_1.fq.gz | wc -l) / 4 ))
#zcat: to decompress .gz file
#wc-1: counts the number of lines in decompressed output
#divsion by 4 to calculate total number of reads (since each fastq entry will have 4 lines)
#echo to show the result of the calculation

#or

#1)
less LGS1_sub_1.fq.gz | wc -l | awk '{print $1 / 4}'

#2) Reads shorter than 150 BP
zcat BEN_CI16_sub_1.fq.gz | awk 'NR % 4 == 2 {print length($0)}' | awk '$1<150' | wc -l
## NR % 4 == 2: This will check every 2nd line of the each four line of FASTAQ and length($0): This will returns the length of the sequence.

#3) Reads longer than 150 BP
zcat BEN_CI16_sub_1.fq.gz | awk 'NR % 4 == 2 {print length($0)}' | awk '$1>150' | wc -l 

#4) Reads not equal to 150 BP
zcat file.fastq.gz | awk 'NR % 4 == 2 {if (length($0) != 150) print length($0)}' | wc -l
#if (length($0) != 150): Checks if the length of the sequence is not equal to 150.
#print length($0): Prints the length of sequences that are not 150.

#5) Reads contaminated with illumination adapters (Adapter sequence: CTGTCTCTTATACACATCT)
#Illumina sequencing can read into adapter sequences if the insert is shorter than the read length, so sequencing reads past the DNA fragment into the adapter. Or due to improper adapter trimming.





