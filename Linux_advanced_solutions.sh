#### Download the datasets from https://zenodo.org/records/14258052

### Structure of fastq files
# fastq files are the sequencing output files files obtained after demultimplexing data from illumina seqeuncers. Most sequencing companies provide this file format for illumina data.
# It consists of four lines


@SRR15369215.126490887 #sequence identifier 
GGACCTTCTGTCATTTCACTCCTTCTGAAGTAAGGAGTGAAGTAAACACGAAGTAAACACGACAGGTTAGTCCTATTCCTTCAAGCAGGAGTACAGAAAAGAATGCAAATTCTGGGTTCTAGCCCAGCTTTTACTCCTATGGTTCTATTT #sequence
+ #seperator
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFAFJJJJJJJJJJJJJJJJFJJJJJJJJJFJJJJJJJJJFFFJJJJJJJJJJJ #base call quality scores (ASCII characters)

#After 100-150 base pairs, the accuracy of base calling decreases. Trimming is usually done to retain only high-quality bases.

#Note: after downloading all the .gz files, they might get saved something as 'BEN_CI16_sub_1.fq.gz?download=1'. To appropriately rename all the .gz files in a single-go, use the following command:
 for file in *.fq.gz?download=1; do
>   mv "$file" "${file%\?download=1}"
> done
#for file in *.fq.gz?download=1: Loops through all files that match the download pattern.
#mv "$file" "${file%\?download=1}": Renames each file by stripping the ?download=1 using ${file%\?download=1}, which is shell syntax for pattern trimming.

#finally you will have files named something as BEN_CI16_sub_1.fq.gz

# Solutions

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
#$() runs a command and gets its output. $(( )) does math with that output

#or

#1)
less LGS1_sub_1.fq.gz | wc -l | awk '{print $1 / 4}'

#for all the files,

for file in *.fq.gz; do
  read_count=$(( $(zcat "$file" | wc -l) / 4 ))
  echo "$(basename "$file"): $read_count reads"
done

#*.fq.gz is the wildcard entry for all the files ending with ".fq.gz"
#zcat "$file" → Decompresses the .fq.gz file.
#wc -l → Counts total lines.
#/ 4 → Since 1 read = 4 lines in a FASTQ file.
#echo ... → Prints: filename: <read count> reads

#every file should have 1000000 reads

#2) Reads shorter than 150 BP
zcat BEN_CI16_sub_1.fq.gz | awk 'NR % 4 == 2 {print length($0)}' | awk '$1<150' | wc -l
#awk 'NR % 4 == 2 {print length($0)}': FASTQ format has reads on every 2nd line of each 4-line block. NR % 4 == 2 selects those lines. length($0) gives the length of the read sequence.
#awk '$1<150': Filters only those reads shorter than 150 base pairs.
## NR % 4 == 2: This will check every 2nd line of the each four line of FASTAQ and length($0): This will returns the length of the sequence.

#for all the files,

for file in *.fq.gz; do
   reads=$(zcat "$file" | awk 'NR % 4 == 2 {print length($0)}' | awk '$1<150' | wc -l)
   echo "$(basename "$file"): $reads reads shorter than 150bp"
 done

 #Every file should have 0 reads shorter than 150bp

#3) Reads longer than 150 BP
zcat BEN_CI16_sub_1.fq.gz | awk 'NR % 4 == 2 {print length($0)}' | awk '$1>150' | wc -l 

#for all the files,

for file in *.fq.gz; do
   reads=$(zcat "$file" | awk 'NR % 4 == 2 {print length($0)}' | awk '$1>150' | wc -l)
   echo "$(basename "$file"): $reads reads longer than 150bp"
 done

#BEN_CI16_sub_1.fq.gz: 0 reads longer than 150bp
#BEN_CI16_sub_2.fq.gz: 0 reads longer than 150bp
#BEN_CI18_sub_1.fq.gz: 0 reads longer than 150bp
#BEN_CI18_sub_2.fq.gz: 0 reads longer than 150bp
#BEN_NW10_sub_1.fq.gz: 454124 reads longer than 150bp
#BEN_NW10_sub_2.fq.gz: 454124 reads longer than 150bp
#BEN_NW12_sub_1.fq.gz: 0 reads longer than 150bp
#BEN_NW12_sub_2.fq.gz: 0 reads longer than 150bp
#BEN_NW13_sub_1.fq.gz: 0 reads longer than 150bp
#BEN_NW13_sub_2.fq.gz: 0 reads longer than 150bp
#BEN_SI18_sub_1.fq.gz: 1000000 reads longer than 150bp
#BEN_SI18_sub_2.fq.gz: 1000000 reads longer than 150bp
#BEN_SI19_sub_1.fq.gz: 1000000 reads longer than 150bp
#BEN_SI19_sub_2.fq.gz: 1000000 reads longer than 150bp
#BEN_SI9_sub_1.fq.gz: 0 reads longer than 150bp
#BEN_SI9_sub_2.fq.gz: 0 reads longer than 150bp
#LGS1_sub_1.fq.gz: 0 reads longer than 150bp
#LGS1_sub_2.fq.gz: 0 reads longer than 150bp

#4) Reads not equal to 150 BP
zcat BEN_CI16_sub_1.fq.gz | awk 'NR % 4 == 2 {if (length($0) != 150) print length($0)}' | wc -l
#if (length($0) != 150): Checks if the length of the sequence is not equal to 150.
#print length($0): Prints the length of sequences that are not 150.

#for all the files,
for file in *.fq.gz; do
  count=$(zcat "$file" | awk 'NR % 4 == 2 {if (length($0) != 150) print length($0)}' | wc -l)
  echo "$(basename "$file"): $count reads not equal to 150bp"
done

#BEN_CI16_sub_1.fq.gz: 0 reads not equal to 150bp
#BEN_CI16_sub_2.fq.gz: 0 reads not equal to 150bp
#BEN_CI18_sub_1.fq.gz: 0 reads not equal to 150bp
#BEN_CI18_sub_2.fq.gz: 0 reads not equal to 150bp
#BEN_NW10_sub_1.fq.gz: 454124 reads not equal to 150bp
#BEN_NW10_sub_2.fq.gz: 454124 reads not equal to 150bp
#BEN_NW12_sub_1.fq.gz: 0 reads not equal to 150bp
#BEN_NW12_sub_2.fq.gz: 0 reads not equal to 150bp
#BEN_NW13_sub_1.fq.gz: 0 reads not equal to 150bp
#BEN_NW13_sub_2.fq.gz: 0 reads not equal to 150bp
#BEN_SI18_sub_1.fq.gz: 1000000 reads not equal to 150bp
#BEN_SI18_sub_2.fq.gz: 1000000 reads not equal to 150bp
#BEN_SI19_sub_1.fq.gz: 1000000 reads not equal to 150bp
#BEN_SI19_sub_2.fq.gz: 1000000 reads not equal to 150bp
#BEN_SI9_sub_1.fq.gz: 0 reads not equal to 150bp
#BEN_SI9_sub_2.fq.gz: 0 reads not equal to 150bp
#LGS1_sub_1.fq.gz: 0 reads not equal to 150bp
#LGS1_sub_2.fq.gz: 0 reads not equal to 150bp

#5) Reads contaminated with illumination adapters (Adapter sequence: CTGTCTCTTATACACATCT)
#Illumina sequencing can read into adapter sequences if the insert is shorter than the read length, so sequencing reads past the DNA fragment into the adapter. Or due to improper adapter trimming.

zcat BEN_CI16_sub_1.fq.gz | awk 'NR % 4 == 2' | grep -c 'CTGTCTCTTATACACATCT'
#zcat BEN_CI16_sub_1.fq.gz: uncompresses the file.
#awk 'NR % 4 == 2': picks only the sequence lines from each read (line 2 out of every 4 in a FASTQ file).
#grep -c 'CTGTCTCTTATACACATCT': counts how many sequences contain the adapter.

#for all the files:

for file in *.fq.gz; do
>     count=$(zcat "$file" | awk 'NR % 4 == 2' | grep -c 'CTGTCTCTTATACACATCT')
>     echo "$file: $count"
> done

#All the files will have 0 contaminations.

