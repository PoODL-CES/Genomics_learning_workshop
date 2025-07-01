
#### the filtered vcf generated after variant filtering may be used for this exercise - 
# machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB.vcf
     

####PCA

# Firstly, .bed, .bi,m and .fam files need to be generated. from a vcf file. 
#.bed: binary file that stores genotype calls for every individual at very SNP.
#.bim: tab-delimited text file with information about each SNP.
#.fam: tab-delimited text file; contains information about each individual/sample.




conda activate plink
plink --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB.vcf \
  --make-bed \
  --double-id \
  --allow-extra-chr \
  --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB
#plink: converts VCF file to PLINK binary format (.bed, .bim, .fam)
#--vcf <file>: specifies the input VCF file containing the genotype data
#--make-bed: converts the input into PLINK binary format.
#--double-id: Treats the entire VCF sample ID as both Family ID (FID) and Individual ID (IID), which is necessary when IDs include underscores (_) that would otherwise be misinterpreted.
#--allow-extra-chr: allows non-standard chromosome names which are not usually allowed in strict plink parsing
#--out output_file: sets the prefix for all the output files.
# bed bim and bam file would be created with the specified prefix.

#### We have to edit the fam file in order to add the region information. The second column of the fam file (IID) is edited to add information about region. We look at the names in the first column and add region information based on that
#   For examples if the sample name has CI - it is from central india, SI means south india and so on.
# This step is done in order to get proper region information while doing PCA.

awk '{
    id = $1

    if (id ~ /_CI/)       region = "CenIndia"
    else if (id ~ /_SI/)  region = "SouIndia"
    else if (id ~ /_NE/)  region = "NorEasIndia"
    else if (id ~ /_NW/)  region = "NorWesIndia"
    else if (id ~ /_NOR/) region = "NorIndia"
    else if (id ~ /_SU/)  region = "Sunderban"
    else if (id ~ /^LGS/) region = "CenIndia"
    else                  region = "Unknown"

    $2 = region
    print
}' machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB.fam > tmp.fam \
&& mv tmp.fam machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB.fam


#### Edit the bim file to replace the chromosome name E2 with any integer value. We need to do this because admixture does not support non-integer strings.

awk '{$1 = 1; print}' new_test.bim > tmp && mv tmp new_test.bim # replacing E2 with 1 





################# For PCA
conda activate plink

  plink --bfile machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB \
  --pca 5 \
  --allow-extra-chr \
  --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB

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
fam <- read.table("machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB.fam", header = FALSE)
eigenvec <- read.table("machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB", header = FALSE)
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
admixture machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB.bed $K
done

## For example if K=3, admixture assumes there are 3 ancestral populations and estimates how much of each individual's genome comes from each of the 3.
## For every value of K, output files with extensions .K.P and .K.Q will be generated

## We will use the .3.Q file for plotting below.



### to plot the admixture results

conda activate R_env

R #(It will open R prompt)

# Load libraries 

library(ggplot2) 

library(reshape2) 

library(dplyr)  

# Read ADMIXTURE Q file and sample IDs from .fam 

q3 <- read.table("machali_...3.Q")  # Replace with actual full filename 

fam <- read.table("machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB.fam")

sample_ids <- fam$V1 

# Add sample IDs 

q3$ID <- sample_ids 


# OPTIONAL: Add groups/populations (e.g., from file or by pattern) 

# If you have a file with sample ID and population: 

# groups <- read.table("group_info.txt", header = TRUE) 

# q3 <- merge(q3, groups, by.x = "ID", by.y = "SampleID") 

# OR, extract group info from sample ID 

q3$Group <- gsub(".*_(.*)", "\\1", q3$ID)  # example: get last part after underscore 

# Reshape data for ggplot 

q3_long <- melt(q3, id.vars = c("ID", "Group")) 

# Sort individuals by group 

q3_long <- q3_long %>% 

  arrange(Group, ID) %>% 

  mutate(ID = factor(ID, levels = unique(ID))) 

# Plot 

p <- ggplot(q3_long, aes(x = ID, y = value, fill = variable)) + 

  geom_bar(stat = "identity", width = 1) + 

  facet_grid(~Group, scales = "free_x", space = "free_x") + 

  theme_minimal() + 

  labs(x = "Individuals", y = "Ancestry Proportion", title = "ADMIXTURE Plot (K=3)") + 

  theme(axis.text.x = element_blank(), 

        axis.ticks.x = element_blank(), 

        panel.spacing = unit(0.5, "lines"), 

        strip.text.x = element_text(angle = 0, face = "bold"), 

        legend.position = "right") + 

  scale_fill_brewer(palette = "Set1") 
# Save as PNG 

ggsave("admixture_K3_grouped.png", p, width = 12, height = 6, dpi = 300)

################# For heterozygosity

## We will use RTG-tools for estimating heterozygosity per sample (https://www.realtimegenomics.com/products/rtg-tools)
## To install RTG tools

conda create -n rtg-tools -c bioconda rtg-tools

conda activate rtg-tools

rtg vcfstats machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_rmvIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_imiss_0.6_miss_0.6_mid95percentile_noZSB.vcf > machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.vcfstats

conda deactivate
