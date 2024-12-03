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
#then, log into the CES cluster 
fastqc BEN_CI16_sub_2.fq.gz
xdg-open BEN_CI16_sub_2_fastqc.html
exit
#then, log into the CES cluster 
scp uid@10.58.21.46:/home/pranav/BEN_CI16_sub_2_fastqc.html .
#open the file from your home directory 
