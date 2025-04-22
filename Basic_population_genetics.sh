#### Download the datasets from https://zenodo.org/records/15259681
### the filtered vcf generated after variant filtering maybe used for this exercise

#PCA

#firstly, .bed, .bim and .fam files need to be generated. from a vcf file. 
#.bed: binary file that stores genotype calls for every individual at very SNP.
#.bim: tab-delimited text file with information about each SNP.
#.fam: tab-delimited text file; contains information about each individual/sample.

# We will do the convertion in two steps. First we will create a ped file using vcftools --plink option

conda activate vcftools

vcftools --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.recode.vcf \
  --plink \
  --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB

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
#--pca 10: calculates the top 10 principal components
#--allow-extra-chr: allows non-standard chromosome names
#out output_file_pca: sets the prefix for output files (example: output_file_pca.eigenval)
#eigenvalues and eigenvectors would be created
conda deactivate


# To plot the PCA we will use R (https://www.r-project.org)
# To install R
conda create -n r_env -c conda-forge r-base r-ggplot2
# Activate the R environment
conda activate r_env
# Activate R
R
eigenvec_data <- read.table("output_file_pca.eigenvec", header=FALSE)
colnames(eigenvec_data) <- c("FID", "IID", paste("PC", 1:10, sep=""))
head(eigenvec_data)
ggplot(eigenvec_data, aes(x=PC1, y=PC2)) +
  geom_point() +
  labs(x="Principal Component 1", y="Principal Component 2", title="PCA Plot: PC1 vs PC2") +
  theme_minimal()
ggsave("pca_plot.png")
q()
conda deactivate

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


# For running the tutorial on CES server: We exit the server and download the file to our computers to view it. For this we exit the server using the 'exit' command and then run 'scp username@IP_address:~/"path to the file on the remote cluster"/Rplots.pdf .'

################# For ADMIXTURE

conda create -n admixture -c bioconda admixture

for K in {2..4}
do 
admixture machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe0.05_mm0.6_meanDPmid95percentile_imiss0.6_noZoo.bed ${K}
done

# We edit the .fam file with the location names for the samples and then paste the location and sample names to the .Q files created by ADMIXTURE
# The .Q file has the ancestry proportion for each individual from each ancestral population
# We paste the first two columns of the .fam file to the .Q file

awk '{print $1"\t"$2}' > geographical_location.txt
paste geographical_location.txt machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe0.05_mm0.6_meanDPmid95percentile_imiss0.6_noZoo.2.Q > machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe0.05_mm0.6_meanDPmid95percentile_imiss0.6_noZoo.2.K


