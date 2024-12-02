#### Download the datasets from 

### Count the number of sequencing reads in the fastq files
solution: for LGS1_sub_1.fq
grep '^@' LGS1_sub_1.fq | wc -l ### grep '^@' Finds all lines starting with @ in the file.
^@: Matches lines where @ appears at the beginning.
wc -l: Counts the number of matching lines.


