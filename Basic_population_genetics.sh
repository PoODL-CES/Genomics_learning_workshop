#### Download the datasets from https://zenodo.org/records/15263700
### the filtered vcf generated after variant filtering maybe used for this exercise

#PCA

#firstly, .bed, .bim and .fam files need to be generated. from a vcf file. 
#.bed: binary file that stores genotype calls for every individual at very SNP.
#.bim: tab-delimited text file with information about each SNP.
#.fam: tab-delimited text file; contains information about each individual/sample.

# We will do the convertion in two steps. First we will create a ped file using vcftools --plink option

conda activate vcftools

vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB.recode.vcf \
  --plink \
  --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB
conda deactivate

# Then we will create the bed file by inputting the ped file to plink

conda activate plink
plink --file machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB \
  --make-bed \
  --double-id \
  --allow-extra-chr \
  --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB

#plink: converts VCF file to PLINK binary format (.bed, .bim, .fam)
#--vcf <file>: specifies the input VCF file containing the genotype data
#--make-bed: converts the input into PLINK binary format.
#--double-id: Treats the entire VCF sample ID as both Family ID (FID) and Individual ID (IID), which is necessary when IDs include underscores (_) that would otherwise be misinterpreted.
#--allow-extra-chr: allows non-standard chromosome names which are not usually allowed in strict plink parsing
#--out output_file: sets the prefix for all the output files.
#output_file.bed, output_file.bim, and output_file.fam would be created
conda deactivate

################# For PCA
conda activate plink

  plink --bfile machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB \
  --pca 5 \
  --allow-extra-chr \
  --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB

#--bfile output_file: loads the plink binary dataset previously created using --make-bed
#--pca 5: calculates the top 5 principal components
#--allow-extra-chr: allows non-standard chromosome names
#out output_file_pca: sets the prefix for output files (example: output_file_pca.eigenval)
#eigenvalues and eigenvectors would be created
conda deactivate


# To plot the PCA we will use R (https://www.r-project.org)
# Activate the R environment
conda activate R_env

# Activate R
# 
R
fam <- read.table("machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.fam", header = FALSE)
eigenvec <- read.table("machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.eigenvec", header = FALSE)
colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:5))
colnames(fam)[2] <- "Region"
eigenvec$Region <- fam$V2
library(ggplot2)
ggplot(eigenvec, aes(x = PC1, y = PC2, color = IID)) +
geom_point(size = 3, alpha = 0.8) +
theme_minimal() +
labs(title = "PCA Plot by Region", x = "PC1", y = "PC2") 

ggsave("pca_by_region.png")
ggsave("pca_by_region.pdf")

#R: Launches R
#fam <- read.table("output_file.fam", header=FALSE): reads the edited .fam file into the dataframe  
#eigenvec_data <- read.table("output_file_pca.eigenvec", header=FALSE): reads the .eigenvec file into the dataframe
#colnames(eigenvec_data) <- c("FID", "IID", paste("PC", 1:5)): assigns column names; first 2 columns are FID and IID while next 5 columns are PC1 to PC5
#colnames(fam)[2] <- "Region": Renames the second column of the .fam file to "Region"
#eigenvec$Region <- fam$V2: Adds a new "Region" column to the eigenvec data frame using the values from the .fam file’s 2nd column
#head(eigenvec_data): displays the first few rows for the purpose of confirmation
#library(ggplot2): loads ggplot2 
#ggplot(eigenvec, aes(x = PC1, y = PC2, color = IID)) +: tells ggplot to use PC1 on the x-axis and PC2 on the y-axis, and color points by IID.
#geom_point(size = 3, alpha = 0.8) +:Adds points to the plot, size 3, and 80% opaque.
#theme_minimal() +: A clean and simple background theme.
#labs(title = "PCA Plot by Region", x = "PC1", y = "PC2") +:Labels the title and axes.
#scle_color_brewer(palette = "Set2"): Applies a color palette ("Set2") from ColorBrewer
#ggsave("pca_by_region.png"): saves the plot in .png format
#ggsave("pca_by_region.pdf"): saves the plot in .pdf format
#q(): Exits the R console

scp 
#scp: secure copy protocol; copies files betweeen a local and a remote computer.
# .: destination on the local machine. copies to the current directory.
#pca_plot.png and Rplots.pdf would be saved in your local home directory


# For running the tutorial on CES server: We exit the server and download the file to our computers to view it. For this we exit the server using the 'exit' command and then run 'scp username@IP_address:~/"path to the file on the remote cluster"/Rplots.pdf .'

################# For ADMIXTURE

conda create -n admixture -c bioconda admixture
conda activate admixture

for K in {2..4}
do 
admixture machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.bed ${K}
done

# We use the edited .fam file with the location names for the samples and then paste the location and sample names to the .Q files created by ADMIXTURE
# The .Q file has the ancestry proportion for each individual from each ancestral population
# We paste the first two columns of the .fam file to the .Q file

for K in {2..4}
do
paste geographical_location.txt machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.${K}.Q > machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.${K}.xls
done

## to plot the admixture results

################# For heterozygosity

## We will use RTG-tools for estimating heterozygosity per sample (https://www.realtimegenomics.com/products/rtg-tools)
## To install RTG tools

conda create -n rtg-tools -c bioconda rtg-tools

conda activate rtg-tools

rtg vcfstats machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.recode.vcf.gz > machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.vcfstats

conda deactivate
