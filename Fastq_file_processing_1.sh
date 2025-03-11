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

#### Trim the raw fastq reads to remove adapter sequences and low quality bases form the ends
