#### Map trimmed reads using bwa mem
## The datasets from https://zenodo.org/records/14258052 contain the reference genome file and its index files

## The reads that were trimmed to remove adapters and low quality bases need to be mapped to a reference genome
# TASK 1: map the reads to the GCA_021130815.1_PanTigT.MC.v3_genomic.fna reference genome that was downloaded

## The mapping creates a .sam file. This occupies lots of storage space on the server so we should convert it to a .bam file
## .sam is abbreviation of Sequence Alignment Map and .bam is Binary Alignment Map

# TASK 2: Convert sam file to bam 

## The reads are arranged randomly in the sam and bam file. They need to sorted as per the genomic corrdinates they occupy.
## Sorting the reads improves the efficiency of downstream analysis like variant calling etc. by giving access to specific regions of the genome.

# TASK 3: sort the bam files

## PCR replicates and optical replicates need to be removed from analysis as they may falsely infale the depth.
#PCR duplicates (a.k.a PCR duplicates): They are identical reads or read pairs resulting from the overamplification during PCR step of library preparation, and not from the actual biological sample. During sequencing, they appear as though the same fragment was read multiple times.
#some fragments get amplified more randomly or due to better "amplifiability". Poor-quality enzymes or unoptimized PCR conditions may lead to non-specific amplification.
#Optical duplicates: Specifically arisen from sequencing machine artifacts. When a cluster is wrongly interpreted as multiple nearby clusters (due to how it fluoresces), you get multiple reads that look like PCR duplicates but were actually created during imaging.
## This step is to be avoided for ddRAD, RADseq and amplicon sequencing since most reads are PCR duplicates or match other reads at the initial few bases

# TASK 4: Remove or mark the duplicate reads

### INDEX AFTER MARKDUPLICATES
## Indexing allows quick access to specific genomic regions and improve performance of downstream analysis tools

# TASK 5: index the sort mark duplicatd bam file for SNP calling

### Estimate some statistics from the bam files

# TASK: estimate the sequencing depth and percent mapping of the reads to the reference genome

# Look for a section Chromosome-wise coverage” or similar
