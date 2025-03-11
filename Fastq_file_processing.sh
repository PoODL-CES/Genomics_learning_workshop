#### Install conda (https://docs.anaconda.com/miniconda/install/)
### For linux systems

mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh

### After installing, close and reopen your terminal application or refresh it by running the following command:
source ~/miniconda3/bin/activate

### To initialize conda on all available shells, run the following command:
miniconda3/bin/./conda init --all


#### For getting summary of fastq files, install fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Install fastqc
conda create -n fastqc -c bioconda fastqc

## Activate the conda environment
conda activate fastqc

## run fastqc
fastqc *.fq.gz

## Deactivate the conda environment
conda deactivate fastqc

#### Trimming Raw fastq files to remove adapter sequences and low quality bases form the ends
#### Install trim_galore
conda create -n trim-galore -c bioconda trim-galore

## activate the conda environment
conda activate trim-galore

## for DNA seqeunced on illumina sequencers with illumina adapters the following maybe done.
## for other library preparations with different adapters, adapter sequences need to be specified.
## for BGI seq different adapters need to be declared.

trim_galore --paired --illumina BEN_CI16_sub_1.fq.gz BEN_CI16_sub_2.fq.gz

## triming for single fasta file
trim_galore --fastqc BEN_CI16_sub_1.fq.gz

## To chcek the output 
less BEN_CI16_sub_1_trimmed.fq.gz

## For bulk data trimming ##
for file in BEN_CI18_sub_1.fq.gz BEN_NW13_sub_1.fq.gz BEN_SI9_sub_1.fq.gz BEN_NW10_sub_1.fq.gz BEN_SI18_sub_1.fq.gz LGS1_sub_1.fq.gz BEN_NW12_sub_1.fq.gz BEN_SI19_sub_1.fq.gz; do
    paired_file=${file/_1.fq.gz/_2.fq.gz}
    trim_galore --paired "$file" "$paired_file"
done
file List:
The loop iterates only over _1.fq.gz files (read 1).
Example: BEN_CI18_sub_1.fq.gz, BEN_NW13_sub_1.fq.gz, etc.
paired_file Variable:

Replaces _1.fq.gz with _2.fq.gz to find the corresponding pair.
trim_galore Command:

Trims both read 1 ($file) and read 2 ($paired_file) together.

#2
for file in *_1.fq.gz; do
  trim_galore --paired "$file" "${file/_1.fq.gz/_2.fq.gz}"
done
# "$file" refers to the current forward read (read 1).
# ${file/_1.fq.gz/_2.fq.gz} dynamically generates the corresponding reverse read (read 2) by replacing _1 with _2 in the file name.

## Deactivate the conda environment
conda deactivate trim-galore
