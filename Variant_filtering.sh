#### Download the vcf file to be used from this link https://zenodo.org/records/15173226


## Applying base quality, genotype quality and hwe filters
## Remove Indels
## Apply the individual missingness filter
## Remove individuals which have missing proportion greater than 60 percent.
## Apply a missingness filter to only keep variants which are common among a proportion of individuals.
## Extract the number of variants remaining upon applying missingness filter from 0.1 to 0.9 and make a .txt file.
## Plot the data in the variant_counts.txt using ggplot2
## Apply max missing 0.6 filter. the one with atleast 60% of data have been retained 
## Use site mean depth filter to get informaion about mean depth at each variant site
## find the percentile depth between 0.05 and 0.975, the middle 95 percentile.
## Apply filters to keep only the variants that have depth values lying in the above range
## perform PCA
