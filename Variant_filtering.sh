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
    count=$(grep -vc "^#" filtered_max_missing_${miss}.recode.vcf)
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

conda deactivate


##################################################################


#perform PCA

#firstly, .bed, .bim and .fam files need to be generated. from a vcf file. 
#.bed: binary file that stores genotype calls for every individual at very SNP.
#.bim: tab-delimited text file with information about each SNP.
#.fam: tab-delimited text file; contains information about each individual/sample.

conda create -n plink -c bioconda plink
conda activate plink
plink --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile.recode.vcf \
  --make-bed \
  --double-id \
  --allow-extra-chr \
  --out output_file

#plink: converts VCF file to PLINK binary format (.bed, .bim, .fam)
#--vcf <file>: specifies the input VCF file containing the genotype data
#--make-bed: converts the input into PLINK binary format.
#--double-id: Treats the entire VCF sample ID as both Family ID (FID) and Individual ID (IID), which is necessary when IDs include underscores (_) that would otherwise be misinterpreted.
#--allow-extra-chr: allows non-standard chromosome names which are not usually allowed in strict plink parsing
#--out output_file: sets the prefix for all the output files.
#output_file.bed, output_file.bim, and output_file.fam would be created

  plink --bfile output_file \
  --pca 10 \
  --allow-extra-chr \
  --out output_file_pca

#--bfile output_file: loads the plink binary dataset previously created using --make-bed
#--pca 10: calculates the top 10 principal components
#--allow-extra-chr: allows non-standard chromosome names
#out output_file_pca: sets the prefix for output files (example: output_file_pca.eigenval)
#eigenvalues and eigenvectors would be created

conda create -n r_env r_base
conda activate r_env
R
install.packages("ggplot2")
library(ggplot2)      
eigenvec_data <- read.table("output_file_pca.eigenvec", header=FALSE)
colnames(eigenvec_data) <- c("FID", "IID", paste("PC", 1:10, sep=""))
head(eigenvec_data)
ggplot(eigenvec_data, aes(x=PC1, y=PC2)) +
  geom_point() +
  labs(x="Principal Component 1", y="Principal Component 2", title="PCA Plot: PC1 vs PC2") +
  theme_minimal()
ggsave("pca_plot.png")
q()

scp username@IP_address:~/"path to the file on the remote cluster"/Rplots.pdf .
#R: Launches R
#install.packages("ggplot2") library(ggplot2): installs and loads ggplot2 for plotting
#eigenvec_data <- read.table("output_file_pca.eigenvec", header=FALSE): reads the .eigenvec file into the dataframe
#colnames(eigenvec_data) <- c("FID", "IID", paste("PC", 1:10, sep="")): assigns column names; first 2 columns are FID and IID while next 10 columns are PC1 to PC10
#head(eigenvec_data): displays the first few rows for the purpose of confirmation
#ggplot(eigenvec_data, aes(x=PC1, y=PC2)) + geom_point() + labs(x="Principal Component 1", y="Principal Component 2", title="PCA Plot: PC1 vs PC2") + theme_minimal()
       # Actually builds the graph
#ggsave("pca_plot.png"): saves the plot in .png format
#q(): Exits the R console

#scp: secure copy protocol; copies files betweeen a local and a remote computer.
# .: destination on the local machine. copies to the current directory.
#pca_plot.png and Rplots.pdf would be saved in your local home directory
