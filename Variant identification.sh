### Variant calling using strelka
conda create -n strelka -c bioconda strelka
#creates a new Conda environment named strelka and installs Strelka2 from the bioconda channel.

conda activate strelka
#Activates strelka

#Index the reference genome file (skip it if you already have .fai index file)
samtools faidx GCA_021130815.1_PanTigT.MC.v3_genomic.fna
#faidx stands for FASTA indexing.

#Strelka is run in a 2 step procedure

# Step:1 - Configuration - to specify the input data and any options pertaining to the variant calling methods themselves

which configureStrelkaGermlineWorkflow.py
#If Strelka2 is installed properly, this command should return the full path

/home/your_username/miniconda3/envs/strelka/bin/configureStrelkaGermlineWorkflow.py --bam name_of_the_file_aligned_reads_sorted_deduplicated.bam --referenceFasta reference_genome_filename.fna --runDir strelka_germline
#/home/your_username/miniconda3/envs/strelka/bin/configureStrelkaGermlineWorkflow.py: This is the Strelka configuration script used for germline variant calling. It prepares the pipeline and generates the necessary files to run the variant calling workflow.
# --bam name_of_the_file_aligned_reads_sorted_deduplicated.bam: Specifies the input BAM file (aligned, sorted, and deduplicated reads). This is the sequencing data from which germline variants will be called.
#  --referenceFasta reference_genome_filename.fna Provides the reference genome in FASTA format. This is the genome against which the BAM file was aligned.
#--runDir strelka_germline: Specifies the output directory where the workflow will be set up. Inside this directory, Strelka will create scripts and configuration files needed to run the variant calling pipeline.

# Step:2 - Workflow Execution - to specify parameters pertaining to how strelka is executed

strelka_germline/runWorkflow.py -m local -j 8
#strelka_germline/runWorkflow.py: This script was generated during the configureStrelkaGermlineWorkflow.py step. It is the main script that runs the variant calling process.
#-m local: Specifies that the workflow will be executed on the local machine.
#-j 8: Defines the number of parallel threads (CPUs) to use. 

# For loop command for multiple files: takes in bam files iteratively and performs configuration and execution steps inside the loop

for bam in *_aligned_reads_sorted_deduplicated.bam; do  # loops through the required bam files
    sample_name=$(basename "$bam" _aligned_reads_sorted_deduplicated.bam) # extracts the sample name from the .bam filename
    output_dir="strelka_germline_${sample_name}" 
    
    mkdir -p "$output_dir"  # make a directory with the sample name extracted in the above steps

    /home/nithinka/miniconda3/envs/strelka/bin/configureStrelkaGermlineWorkflow.py \
        --bam "$bam" \
        --referenceFasta GCA_021130815.1_PanTigT.MC.v3_genomic.fna \
        --runDir "$output_dir"                                                      # run configuration step

    "$output_dir/runWorkflow.py" -m local -j 8                                      # workflow execution step

done

# an output directory will be present for all input bam files containing results of variant calling.

### Variant filtering

#1) filtering passed from non passed
#using bcftools

conda install -c bioconda bcftools
bcftools --version (#this is confirmatory step)
conda activate bioinfo
bcftools view -f PASS -o passed_variants.vcf.gz "input_file_name".vcf.gz

#bcftools: Calls the bcftools program, a widely used tool for processing VCF/BCF files.
#view: Opens and processes the VCF file.
#-f PASS: Filters out variants and keeps only those with "PASS" in the FILTER column.
#-o: passed_variants.vcf.gz	Specifies the output file name where filtered variants will be saved.
#"input_file_name".vcf.gz: The name of the input compressed VCF file containing variant calls.

#2) filtering out indels
#to avoid alignment issues
bcftools view -v snps -o snps_only.vcf.gz passed_variants.vcf.gz

#3) Filtering out Minor allele counts
 # a filter to remove very rare variants, which might be due to sequencing errors rather than real genetic variation. Rare variants are often sequencing errors rather than real mutations.
 bcftools view -i 'MAC >= 3' -o mac_filtered.vcf.gz snps_only.vcf.gz

#4) genotype quality filter
#ensure confidence in genotype calls
bcftools view -i 'FMT/GQ >= 30' -o gq_filtered.vcf.gz mac_filtered.vcf.gz

#5) base quality filter
#removes low quality base calls
bcftools view -i 'QUAL >= 30' -o bq_filtered.vcf.gz gq_filtered.vcf.gz

#6) count the number of SNP's
bcftools view -H -v snps bq_filtered.vcf.gz | wc -l

#### APPLYING VARIANT FILTERS ON A VCF FILE

# Install vcftools in a new conda environment
conda create -n vcftools -c bioconda vcftools
conda activate vcftools

# We have been provided with a vcf file with few filters already applied
# filename - machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz
# when applying new filters, change the name of the output file to reflect the new filters that have been applied.

## Applying base quality, genotype quality and hwe filters:
vcftools --gzvcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz /
--minQ 30 --minGQ 30 --hwe 0.05  --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05 --recode

# --minQ  : Minimum base quality - Filters out variant sites where the supporting base calls have low Phred quality scores, reducing the chance of including sequencing errors.
# --minGQ : Minimum genotype quality - Filters out genotypes with low confidence in the assigned genotype, ensuring that only reliable genotype calls are retained.
# --hwe   : Filter out variants which are under selection / deviate fromm hardy weinberg equillibrium. Removes sites with p values less than the set threshold.

## Apply a missingness filter to only keep variants which are common among a proportion of individuals.

## Do a for loop to apply missingness filter from 0.1 to 0.9 (variants which are found in 10% to 90%)

for p in {10..90..10}; do
  perc=$(echo "scale=2; $p / 100" | bc)
  vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05.recode.vcf \
           --max-missing $perc \
           --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_miss${p} \
           --recode
done
# a different output file will be generated for each value of the missingness filter.

## Extract the number of variants remaining upon applying missingness filter from 0.1 to 0.9 and make a .txt file.

for i in {10..90..10}; do
    count=$(grep -vc "^#" filtered_max_missing_${i}.recode.vcf)
    echo "${i} ${count}" >> variant_counts.txt
done

#### Plot the data in the variant_counts.txt using ggplot2

## create a new conda environment ggplot2
conda create -n r_env r-base
conda activate r_env
R --version
conda create -n ggplot2 -c conda-forge r-ggplot2
conda activate ggplot2

#pranav
R -e "library(ggplot2); library(scales); data <- read.table('variant_counts.txt', header=FALSE, skip=1, col.names=c('Missingness', 'Variants')); data\$Missingness <- as.numeric(as.character(data\$Missingness)); data\$Variants <- as.numeric(as.character(data\$Variants)); p <- ggplot(data, aes(x=Missingness, y=Variants)) + geom_line(color='blue') + geom_point(color='red') + ggtitle('Number of Passed Variants vs. Missingness Filter') + xlab('Max Missingness (%)') + ylab('Number of Variants') + scale_y_continuous(labels = comma) + scale_x_continuous(limits = c(10, 90)) + theme_minimal(); ggsave('variant_plot.png', plot=p)"

# plot the data in the variants.txt file (Nithin)
 R -e "library(ggplot2); library(scales); data <- read.table('variant_counts.txt', header=FALSE, col.names=c('Missingness', 'Variants')); p <- ggplot(data, aes(x=Missingness, y=Variants)) + geom_line(color='blue') + geom_point(color='red') + ggtitle('Number of Passed Variants vs. Missingness Filter') + xlab('Max Missingness (%)') + ylab('Number of Variants') + scale_y_continuous(labels = comma) + scale_x_continuous(limits = c(10, 90)) + theme_minimal(); ggsave('variant_plot.png', plot=p)"

# this will generate a variant_plot.png which can be visualized.



