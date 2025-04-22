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

# For joint variant calling using all samples

    configureStrelkaGermlineWorkflow.py \
        --bam indv1.bam \
        --bam indv2.bam \
        --bam indv3.bam \
        --referenceFasta GCA_021130815.1_PanTigT.MC.v3_genomic.fna \
        --runDir output_folder                                                     # run configuration step

    "$output_dir/runWorkflow.py" -m local -j 8                                      # workflow execution step

done

# an output directory will be present for all input bam files containing results of variant calling.

