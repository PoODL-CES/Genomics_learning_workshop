#### The raw fastq reads sometimes have adapter seqeunces and there are low quality bases towards the ends which needs to be removed
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

#### Obtain summary of fastq files. Count the number of reads in fastq file, the distribution of read lengths and quality
#### For getting summary of fastq files, install fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Install fastqc
conda create -n fastqc -c bioconda fastqc
#conda create: creats a new conda environment. Conda environments are isolated environments that allows to manage dependencies and packages for specific projects without affecting the global setup.
#-n: name of the environment (here it would be fastqc).
#-c: specifices the channel used to download packages (here we are getting fastqc packages in the bioconda channel).
#fastqc: name of the package or the tool we are installing in the new environment. Used for quality control of sequencing data.

## Activate the conda environment
conda activate fastqc

## run fastqc
fastqc *.fq.gz

## Output files generated after running fastqc = "BEN_CI16_sub_1_fastqc.zip" and "BEN_CI16_sub_1_fastqc.html".
# .html file includes visualisations and details about the quality metrics of the sequencing data.
# .gzip file includes summary of the main quality control metrics, detailed quality control data in text format, graphs and images using in .html report.
# .gzip would have much more detailed information

#for simplicity; you can transfer all the .zip files in a global zip folder and all the .html files in a global html folder.

mkdir zip_files html_files #"zip_files" and "html_files" directories would be created
mv *.zip zip_files/ #all the .zip files would be moved to "zip_files"
mv *.html html_files/ #all the .html files would be moved to "html_files"

#to visualize the .html files, they can be copied locally to your local computer from the remote cluster

exit #logout from the remote cluster
scp -r name@IP_address:path to the file or directory ~/

## Deactivate the conda environment
conda deactivate fastqc

#### Trim the raw fastq reads to remove adapter sequences and low quality bases form the ends
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
less BEN_CI16_sub_1_val_1.fq.gz

## For bulk data trimming ##

# OPTION 1: Name all files in the following command
for file in BEN_CI18_sub_1.fq.gz BEN_NW13_sub_1.fq.gz BEN_SI9_sub_1.fq.gz BEN_NW10_sub_1.fq.gz BEN_SI18_sub_1.fq.gz LGS1_sub_1.fq.gz BEN_NW12_sub_1.fq.gz BEN_SI19_sub_1.fq.gz; do
    paired_file=${file/_1.fq.gz/_2.fq.gz}
    trim_galore --paired "$file" "$paired_file"
done
#file List:
#The loop iterates only over _1.fq.gz files (read 1).
#Example: BEN_CI18_sub_1.fq.gz, BEN_NW13_sub_1.fq.gz, etc.
#paired_file Variable:
#Replaces _1.fq.gz with _2.fq.gz to find the corresponding pair.
#trim_galore Command:
#Trims both read 1 ($file) and read 2 ($paired_file) together.

#2 OPTION 2: Use * as wildcard to consider any file with the format *_1.fq.gz as an input
for file in *_1.fq.gz; do
  trim_galore --paired "$file" "${file/_1.fq.gz/_2.fq.gz}"
done
# "$file" refers to the current forward read (read 1).
# ${file/_1.fq.gz/_2.fq.gz} dynamically generates the corresponding reverse read (read 2) by replacing _1 with _2 in the file name.

# Output files generated after running trim-galore : 
#BEN_CI16_sub_1.fq.gz_trimming_report.txt (summary report generated during the trimming process)
#BEN_CI16_sub_1_val_1.fq.gz (val_1: refers to a validated or quality-checked file. In tools like Trim Galore, files are renamed with val_1 post trimming process)

## Deactivate the conda environment
conda deactivate trim-galore
