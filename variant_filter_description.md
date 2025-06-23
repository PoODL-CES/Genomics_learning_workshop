| Filter Name       | Description                                                                                                                                                                            |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `passOnly`        | Retained only variants that passed all internal quality filters from the variant caller (Strelka in this case).                                                                        |
| `biallelicOnly`   | Kept only biallelic variants (i.e., sites with exactly two alleles) for downstream compatibility and clarity.                                                                          |
| `rmvIndels`       | Removed insertions and deletions (INDELs), keeping only SNPs (single nucleotide polymorphisms).                                                                                        |
| `minMAF0Pt05`     | Filtered for variants with a **minor allele frequency (MAF) ≥ 0.05**, ensuring informative polymorphisms and removing rare variants.                                                   |
| `chr_E2`          | Restricted to a specific chromosome/region of interest (e.g., `chr_E2`).                                                                                                               |
| `minDP3`          | Required a **minimum depth (DP) of 3** per genotype to avoid false positives due to low coverage.                                                                                      |
| `minQ30`          | Ensured that each site has a **minimum site quality score of 30**, reflecting high confidence in the variant.                                                                          |
| `minGQ30`         | Filtered genotypes to retain only those with a **Genotype Quality (GQ) ≥ 30**, removing uncertain genotype calls.                                                                      |
| `hwe_0.05`        | Removed variants deviating significantly from **Hardy-Weinberg Equilibrium** (p < 0.05), which could indicate genotyping errors or population substructure.                            |
| `imiss_0.6`       | Filtered out individuals with **>60% missing genotypes** to maintain data quality.                                                                                                     |
| `miss_0.6`        | Removed variants missing in **>60% of individuals**, ensuring sufficient representation across samples.                                                                                |
| `mid95percentile` | Likely indicates removal of extreme outliers in genotype depth or quality, retaining the **middle 95%** of data (please clarify if this is custom logic).                              |
| `noZSB`           | Custom filter|


## Minor Allele Frequency (MAF)

**MAF** is the frequency at which the **less common allele** occurs at a genetic locus in a population.
 *Example:* If allele A appears in 90% of individuals and G in 10%, then **MAF = 0.10**.

---

#### Why apply `minMAF0Pt05`?

This filter removes variants with **MAF < 0.05**, retaining only those where the minor allele is present in **at least 5%** of individuals.

**Benefits**:

* Keeps informative, common variants
* Reduces noise from rare alleles or sequencing errors
* Improves reliability for analyses like PCA, GWAS, or Fst

---


## Depth Filter: minDP3
DP (Depth of Coverage) refers to the number of sequencing reads that support a variant call at a given site in a sample.

#### What does minDP3 mean?
This filter ensures that each genotype (or site) is supported by at least 3 reads. Variants with DP < 3 are excluded.

**Benefits**:

* Ensures reliable variant calls
* Reduces errors from low-coverage or spurious reads
* Improves accuracy in downstream analyses

 ---

## Quality Filter: minQ30
This filter retains only variants with a site quality score ≥ 30, indicating high confidence in the variant call. The site quality score (QUAL) in a VCF file reflects the statistical confidence that a variant is real (not a sequencing error).

It's usually Phred-scaled, meaning:

QUAL 30 = 99.9% confidence the variant is true

Higher scores = higher reliability

**Benefits**:
* Removes low-quality or uncertain variants
* Reduces false positives
* Improves overall reliability of the dataset

---
## Genotype Quality Filter: minGQ30
The genotype quality (GQ) score represents the confidence in the individual genotype call at a variant site. It's Phred-scaled, so:

GQ 30 = 99.9% confidence that the genotype is correct

**Benefits**:
* Removes low-confidence genotype calls
* Reduces false genotyping results
* Improves the accuracy of downstream sample-level analyses

---

## Hardy-Weinberg Equilibrium (HWE)
The HWE test checks whether the observed genotype frequencies match expected frequencies under random mating.

The filter hwe_0.05 removes variants that deviate significantly from equilibrium (p-value < 0.05).

**Why?**

* Variants failing HWE may indicate genotyping errors, population structure, or selection
* Helps ensure biological validity of variants used in population-level studies


---


##  Missingness per Sample: imiss_0.6
This filter removes samples with >60% missing genotype data.

** Why 60%?

* Strikes a balance between data quality and dataset retention
* Removes samples or sites with excessive missing genotypes, which can skew analyses
* Still retains enough data for statistical power in most population-level studies

 ** Why not 70% or 80%?

* Looser thresholds allow too much missing data, which can:
* Distort PCA or clustering
* Reduce confidence in allele frequency estimates
* Introduce noise in GWAS or selection scans

Choosing 60% is a moderate and commonly used cutoff that helps ensure data reliability without being overly restrictive.

**Benefits**

* Highly incomplete samples can bias analyses like PCA or ADMIXTURE
* Improves dataset reliability by keeping well-covered samples


---

##  Missingness per Site: miss_0.6
This removes variant sites where more than 60% of samples have missing data.

**Benefits**
* Ensures you're analyzing variants seen in most samples
* Reduces noise from poorly genotyped regions


---


## What is mid95percentile?
This filter retains only those variants where the depth per site falls within the middle 95% of the overall depth distribution.


Using the middle 95 % of the depth distribution is a standard “trim the tails” approach:

### Reason	what it does

* Removes extreme outliers	Cuts the lowest 2.5 % (often under-covered, error-prone sites) and the highest 2.5 % (often multi-mapped or repetitive regions with artificially high depth).

* Retains the vast majority of data	95 % keeps almost everything that is biologically informative while discarding only the most suspicious extremes.

* Balances stringency and power	Tighter cutoffs (e.g., 90 %) risk losing useful data; looser ones (e.g., 99 %) leave more outliers that can skew statistics such as allele frequencies or depth-based filters.

In short, the 95 % window is a pragmatic compromise: it preserves most variants but shields downstream analyses from depth-related artifacts.

**Benefits**
* Removes extremely low-depth sites that may be unreliable
* Removes very high-depth sites, which could indicate repetitive regions or mapping artifacts
* Focuses the dataset on typical, well-covered variants, improving data quality and consistency

 It effectively trims the top and bottom 2.5% of sites based on depth.

---
