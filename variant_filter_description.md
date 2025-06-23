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
| `noZSB`           | Custom filter — possibly removes variants in **Z-score Signal Bias (ZSB)** regions or those affected by strand bias or alignment artifacts (you may want to describe this explicitly). |


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
