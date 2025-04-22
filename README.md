## This hands-on workshop is designed to introduce students to the basic analysis pipeline for Next-Generation Sequencing (NGS) reads, specifically from the Illumina platform for population genetics analysis.  
We’ll focus on whole genome resequencing data (WGS), though the same pipeline may also be applied to RAD-seq or amplicon-seq data (with minor modifications) if a reference genome is available.

By the end of the workshop, participants will understand how to:

1. Work with raw sequencing data (`FASTQ` files)  
2. Perform quality control and preprocessing  
3. Map reads to a reference genome  
4. Call and filter genetic variants (SNPs/indels)

---

###  Pre-requisites & Setup  
Before joining the workshop, please make sure you can access a Unix/Linux shell environment, as all commands and tools will be run from the command line.  
It is a pre-requisite for the workshop to install some form of Linux/shell terminal. For example, those using Windows may try [MobaXterm](https://mobaxterm.mobatek.net), and Mac/Linux systems should already have a terminal available.  
Go through the file `Understand_Linux_commands` to get a good overview of the command language we will be using in this workshop.


---

###  To get the most out of this workshop, we recommend checking out the following resources:

**How Illumina Sequencing works?**  
Check out this brief and engaging video to understand the sequencing technology behind your data:  
[https://www.youtube.com/watch?v=fCd6B5HRaZ8&t=238s](https://www.youtube.com/watch?v=fCd6B5HRaZ8&t=238s)

---

**Linux Shell basics:**  
Start here if you're new to command-line tools.  

- Basic Tutorial: [`Linux_basics.sh`](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Linux_basics.sh)  
- Solutions File: [`Linux_basics_solutions.sh`](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Linux_basics_solutions.sh)

---

**Advanced Linux Commands:**  
For those who have completed the above tutorial or are already familiar with the basics and want to go deeper to understand further downstream steps.

- Advanced Tutorial: [`Linux_advanced.sh`](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Linux_advanced.sh)  
- Solutions File: [`Linux_advanced_solutions.sh`](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Linux_advanced_solutions.sh)

---

### Dataset Download Instructions:

We’ll be using real WGS resequencing data for this workshop. The dataset is a small subset for demonstration purposes.  
**Download Here**: [Download Here](https://zenodo.org/records/14258052)
 

Instructions to download the data from the website:  
1. Go to the link above.  
2. Make a separate directory to download your files using the command `mkdir name_of_directory` in your Linux shell.  
3. Copy the link and download the dataset folder using the `wget` command on your Linux shell.

---

### Workshop Pipeline Breakdown:

Below is the step-by-step pipeline that we’ll be learning and executing during the workshop:

---

**1. Initial Processing of FASTQ Files**  
Script: [`Fastq_file_processing_solutions.sh`](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Fastq_file_processing_solutions.sh)
  
This step involves:  
a) Checking the quality of raw reads (using tools like `fastqc`)  
b) Trimming adapters or low-quality bases using tools such as `Trim-galore` or `trimmomatic`  
c) Producing cleaned `FASTQ` files ready for mapping  

 *Why it's important:* Poor-quality reads can lead to false variant calls or poor mapping, so quality control is crucial.

---

**2. Mapping Reads to Reference Genome**  
Script: [`Read_mapping_and_sorting_solutions.sh`](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Read_mapping_and_sorting_solutions.sh)
  
Here we:  
a) Use an aligner like `BWA` or `Bowtie2` to map the cleaned reads to a reference genome  
b) Convert output from `SAM` to `BAM` format (compressed and binary)  
c) Sort the `BAM` files by genomic coordinates (using `samtools sort`)  
d) Index the `BAM` files so tools can access them efficiently  

 *Why it's important:* Proper alignment is the foundation for all downstream analysis. Sorting/indexing ensures quick and efficient variant detection.

---

**3. Identifying Variants (SNPs/Indels)**  
Script: [`Variant_identification.sh`](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Variant_identification.sh)
  
In this step:  
a) Call variants (Single Nucleotide Polymorphisms and Indels) using tools like `bcftools` or `strelka`  
b) Produce a `VCF` file (Variant Call Format) with detailed information about each variant  
c) Optionally apply variant filtering to remove false positives or low-confidence variants  

 *Why it's important:* This is the core step of population genetics—identifying the genetic differences across samples.

---

 **4. Filtering Variants**  
Script: [`Variant_filtering.sh`](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Variant_filtering.sh)
  
This step involves:  
a) Applying quality-based filters on the raw VCF file using tools like `bcftools` or `vcftools`  
b) Removing variants with low depth, poor quality scores, or those that don't meet desired thresholds  
c) Producing a filtered, high-confidence VCF file for downstream analysis  

*Why it's important:* Filtering improves the reliability of variant calls and ensures that only biologically relevant, high-confidence variants are retained for further analysis.

---

**5. Population Structure & Visualization**  
Script: [`Basic_population_genetics.sh`](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Basic_population_genetics.sh)

This step includes:  
a) Converting filtered VCF files into formats suitable for statistical analyses (e.g.,`.bed`, `.ped`) using tools like `plink`  
b) Performing Principal Component Analysis (PCA) to visualize population structure and sample clustering  
c) Interpreting population-level genetic structure from the PCA plots  

*Why it's important:* PCA helps in understanding the underlying genetic variation across samples and can highlight population-level patterns or sample outliers.

