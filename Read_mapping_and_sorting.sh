#### Map reads using bwa mem
### Install bwa
conda create -n bwa -c bioconda bwa
conda activate bwa

bwa mem GCA_021130815.1_PanTigT.MC.v3_genomic.fna ../fq_files/BEN_NW10_sub_1_val_1.fq.gz ../fq_files/BEN_NW10_sub_2_val_2.fq.gz > BEN_NW_10_aligned_reads.sam


samtools view -S -b BEN_NW_10_aligned_reads.sam > BEN_NW_10_aligned_reads.bam
samtools sort BEN_NW_10_aligned_reads.bam -o BEN_NW_10_sorted_reads.bam
samtools index BEN_NW_10_sorted_reads.bam
samtools flagstat BEN_NW_10_sorted_reads.bam
