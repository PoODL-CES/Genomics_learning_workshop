# Answer key — Questions 1–4

# 1) Compare PCA clustering with ADMIXTURE (PC1 vs PC2, ADMIXTURE K=4)

**How to produce figures**

* PCA (PLINK):

```bash
plink --bfile prefix --pca 10 --out prefix
# produces prefix.eigenvec and prefix.eigenval
```

* ADMIXTURE:

```bash
admixture --cv prefix.bed 4 | tee log4.txt
# produces prefix.4.Q
```

* Plot using your R scripts (PCA: PC1 vs PC2; ADMIXTURE: prefix.4.Q).

**What to look for**

* **Concordant clusters**: If PCA clusters correspond to ADMIXTURE components, then groups of individuals forming a tight cluster on the PC1–PC2 scatter should also show a single dominant ancestry color in the K=3 barplot.
---

# 2) Create a PC2 vs PC3 plot and interpret structure

**How to produce the plot**

* In R: swap axes or change columns used in ggplot:

```r
p <- ggplot(eigenvec, aes(x = PC2, y = PC3, color = Region)) + geom_point()
ggsave("pca/pca_PC2_PC3.png", p, width=7, height=5)
```
**How to interpret relative to ADMIXTURE**

* If PC2–PC3 exposes a split that corresponds to a subcomponent in ADMIXTURE (e.g., Anc2 splits into two subgroups when K increases), that supports real substructure.
* If PC2–PC3 shows separation not mirrored by ADMIXTURE K=3, try higher K value.
---

# 3) Identify outliers and verify them across both methods

**How to identify outliers in PCA**

* Visual inspection: find points isolated from main clusters.
* Numeric rule-of-thumb: flag samples with |PC1| or |PC2| > 3 standard deviations.
R example to flag:

```r
pc_mean <- colMeans(eigenvec[, c("PC1","PC2")])
pc_sd <- apply(eigenvec[, c("PC1","PC2")], 2, sd)
outliers <- eigenvec[
  abs(eigenvec$PC1 - pc_mean["PC1"]) > 3*pc_sd["PC1"] |
  abs(eigenvec$PC2 - pc_mean["PC2"]) > 3*pc_sd["PC2"], ]
```

**Verify them in ADMIXTURE**

* Check the `prefix.3.Q` rows for those sample IDs:

  * Are they extreme in ancestry fractions? (e.g., ~100% one component) — might be reference/population-specific samples.
  * Are they mixed (~50/50)? — could be recent admixed individuals.
  * Are they unusual compared to their labeled population? — might be mislabeling or migration.

**Re-run PCA without outliers**

1. Remove outlier(s) from the genotype dataset (create new PLINK files) and re-run `plink --pca`.
2. Compare variance explained and cluster tightness.

* If clusters become tighter and PC variance shifts, outliers influenced the original PCA.

---

# 4) Determine optimal K (K = 2–5) and compare to PCA

**How to run ADMIXTURE with CV and record results**

```bash
for K in 2 3 4 5; do
  admixture --cv prefix.bed $K 2>&1 | tee log${K}.txt
done
# Extract CV values
grep -h "CV error" log*.txt
```

* Record CV error for each K. The K with the *lowest CV error* is typically the best-supported model (but consider biological significance too).

**What to report and interpret**

* **Optimal K**: State which K had lowest CV.

  * *Example:* “CV error: K2=0.58, K3=0.41, K4=0.43, K5=0.46 → K=3 lowest → choose K=3.”
* **Compare to PCA**:

  * If PCA shows N well-separated clusters, optimal K should correspond roughly to N.

* **Effect of increasing K**:

  * More K may reveal substructure but can also produce tiny components that are noise or population-specific drift. Watch for components that only appear in 1–2 samples.


**How to justify biological choice of K**

* Use CV as objective metric.
* Inspect barplots across K values: prefer K that yields interpretable, geographically meaningful components rather than many tiny components.
* Cross-check with PCA: number and separation of PCA clusters, and whether additional PCs explain meaningful variation.


---
