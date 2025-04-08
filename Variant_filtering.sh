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

#apply max missing 0.6 filter. the one with atleast 60% of data have been retained

vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05.recode.vcf --max-missing 0.6 --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_mm0.6 --recode

#find the percentile depth between 0.05 and 0.975, the middle 95 percentile.
#below 0.05 because low conf
#above 97.5 because of presence of repetitive regions or less complex regions
#Using R 
conda activate vcftools
conda activate R
#or
#conda activate r_env
R
depth_data <- read.table("out.ldepth.mean", header = TRUE)
depths <- depth_data$MEAN_DEPTH
lower_bound <- quantile(depths, 0.025, na.rm = TRUE)
upper_bound <- quantile(depths, 0.975, na.rm = TRUE)
cat("Middle 95% range:", lower_bound, "to", upper_bound, "\n")
#Middle 95% range: 12.3 to 22.89705

#filter rows with range
mid_95 <- depth_data[depths >= lower_bound & depths <= upper_bound, ]

#save output
write.table(mid_95, file = "mid_95_percentile_loci.txt", row.names = FALSE, quote = FALSE)

#output file:  mid_95_percentile_loci.txt

