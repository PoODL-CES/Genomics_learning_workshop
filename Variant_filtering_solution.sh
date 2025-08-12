### Variant filtering
### Download the vcf file to be used from https://zenodo.org/records/15173226

# Variant filtering can be done using different tools like bcftools, vcftools, gatk etc. Here we will provide the solutions for filtering using vcftools
# Please note that the downloaded vcf file already has a few filters applied



#### APPLYING VARIANT FILTERS ON A VCF FILE

# Install vcftools in a new conda environment
conda create -n vcftools -c bioconda vcftools
conda activate vcftools

# We have been provided with a vcf file with few filters already applied
# filename - machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz
# when applying new filters, change the name of the output file to reflect the new filters that have been applied.

## Applying base quality, genotype quality and hwe filters:
vcftools --gzvcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz \
--minQ 30 --minGQ 30 --hwe 0.05  --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05 --recode

# --minQ  : Minimum base quality - Filters out variant sites where the supporting base calls have low Phred quality scores, reducing the chance of including sequencing errors.
# --minGQ : Minimum genotype quality - Filters out genotypes with low confidence in the assigned genotype, ensuring that only reliable genotype calls are retained.
# --hwe   : Filter out variants which are under selection / deviate fromm hardy weinberg equillibrium. Removes sites with p values less than the set threshold.

### Remove Indels
vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05.recode.vcf --remove-indels \ 
--out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05 --recode

### Apply the individual missingness filter

vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05.recode.vcf --missing-indv

# the above command gives you an output file out.imiss. If you open this file you can see the fraction of missing sites for each inidividual in the column fmiss.
# we want to remove the individuals which have a high proportion of missing sites

## Remove individuals which have missing proportion greater than 60 percent.

 vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05.recode.vcf \
 --remove <(awk '$5 > 0.6 {print $1}' out.imiss) --recode --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6



#### Apply a missingness filter to only keep variants which are common among a proportion of individuals.

## Do a for loop to apply missingness filter from 0.1 to 0.9 (variants which are found in 10% to 90%)

for miss in {10..90..10}; do
  perc=$(echo "scale=2; $miss / 100" | bc)
  vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6.recode.vcf \
           --max-missing $perc \
           --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_${miss} \
           --recode
done
# a different output file will be generated for each value of the missingness filter.

## Extract the number of variants remaining upon applying missingness filter from 0.1 to 0.9 and make a .txt file.

for miss in {10..90..10}; do
    count=$(grep -vc "^#" machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_${miss}.recode.vcf)
    echo "${miss} ${count}" >> variant_counts.txt
done

#### Plot the data in the variant_counts.txt using ggplot2

## create a new conda environment ggplot2
conda create -n ggplot2 -c conda-forge r-ggplot2
conda activate ggplot2

#pranav
R -e "library(ggplot2); library(scales); data <- read.table('variant_counts.txt', header=FALSE, skip=1, col.names=c('Missingness', 'Variants')); data\$Missingness <- as.numeric(as.character(data\$Missingness)); data\$Variants <- as.numeric(as.character(data\$Variants)); p <- ggplot(data, aes(x=Missingness, y=Variants)) + geom_line(color='blue') + geom_point(color='red') + ggtitle('Number of Passed Variants vs. Missingness Filter') + xlab('Max Missingness (%)') + ylab('Number of Variants') + scale_y_continuous(labels = comma) + scale_x_continuous(limits = c(10, 90)) + theme_minimal(); ggsave('variant_plot.png', plot=p)"

# plot the data in the variants.txt file (Nithin)
 R -e "library(ggplot2); library(scales); data <- read.table('variant_counts.txt', header=FALSE, col.names=c('Missingness', 'Variants')); p <- ggplot(data, aes(x=Missingness, y=Variants)) + geom_line(color='blue') + geom_point(color='red') + ggtitle('Number of Passed Variants vs. Missingness Filter') + xlab('Max Missingness (%)') + ylab('Number of Variants') + scale_y_continuous(labels = comma) + scale_x_continuous(limits = c(10, 90)) + theme_minimal(); ggsave('variant_plot.png', plot=p)"

# this will generate a variant_plot.png which can be visualized.

#### Apply max missing 0.6 filter. the one with atleast 60% of data have been retained *

# *If you have performed the above step to apply a for loop for max missing 0.1 to 0.9, skip the below step and use the output file for the 0.6 max missingness filter for subsequent steps.

vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6.recode.vcf --max-missing 0.6 --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6 --recode


#### Use site mean depth filter to get informaion about mean depth at each variant site

vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6.recode.vcf --site-mean-depth --out depth

# Above step will give an output file depth.ldepth.mean


#find the percentile depth between 0.05 and 0.975, the middle 95 percentile.
#below 0.05 because low conf
#above 97.5 because of presence of repetitive regions or less complex regions

conda deactivate 

conda create -n R_env -c conda-forge r-base
conda activate R_env

R
depth_data <- read.table("depth.ldepth.mean", header = TRUE)
depths <- depth_data$MEAN_DEPTH
lower_bound <- quantile(depths, 0.025, na.rm = TRUE)
upper_bound <- quantile(depths, 0.975, na.rm = TRUE)
cat("Middle 95% range:", lower_bound, "to", upper_bound, "\n")
#Middle 95% range: 12.8302 to 24.6604

conda deactivate

### Apply filters to keep only the variants that have depth values lying in the above range

conda activate vcftools

vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6.recode.vcf --min-meanDP 12.8302 --max-meanDP 24.6604 \
--recode --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile

#### Remove zoo samples from the vcf file ( samples starting with id ZSB)

grep 'ZSB' out.imiss | awk '{print $1}'> zsb_samples.txt # extracts sample ids whcih start with ZSB and write it into a text file.

vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile.recode.vcf \
--remove zsb_samples.txt --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB --recode

# take the sample names from zsb_samples.txt and remove those from our vcf file.



##################################################################


