---
layout: post
title: ASCA Analysis Methods for RNAseq and Lipidomics Data
date: '2025-12-19'
categories: Larval_Symbiont_TPC_2023
tags: Mcapitata GeneExpression Lipids Molecular ASCA Methods R
---

This post describes the analytical methods, approach, and main results from three ASCA (ANOVA Simultaneous Component Analysis) analyses conducted on *Montipora capitata* larval multi-omic data: (1) RNAseq gene expression analysis, (2) functional enrichment of RNAseq results, and (3) lipidomics analysis. These analyses examine how temperature and parental phenotype affect molecular phenotypes in coral larvae using a full factorial experimental design.

# Overview

ASCA (ANOVA Simultaneous Component Analysis) is a powerful multivariate analysis approach that combines ANOVA with Principal Component Analysis (PCA) to identify coordinated patterns of variation in high-dimensional datasets. In this study, I applied ASCA to examine molecular responses of *M. capitata* larvae to:

- **Temperature treatments**: 27°C (ambient), 30°C (moderate stress), and 33°C (high stress)
- **Parental phenotypes**: Wildtype (mixed *Cladocopium* and *Durusdinium* symbionts), Bleached (parents with bleaching history; *Cladocopium* only), and Nonbleached (parents without bleaching history; mixed symbionts)

The analyses are implemented in three related scripts:
1. [**rna-seq_ASCA.Rmd**](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_ASCA.Rmd): ASCA analysis of gene expression data
2. [**rna-seq_functional_enrichment_ASCA_topGO.Rmd**](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_functional_enrichment_ASCA_topGO.Rmd): Functional enrichment of genes driving each pattern
3. [**lipids_ASCA.Rmd**](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/lipids_ASCA.Rmd): ASCA analysis of lipidomic data

# Why ASCA is Appropriate for Full Factorial Designs

ASCA is particularly well-suited for analyzing -omics data from full factorial experimental designs for several key reasons:

## Advantages for High-Dimensional Data

1. **Handles more variables than observations**: Unlike traditional MANOVA, ASCA can analyze datasets where the number of variables (genes, lipids) far exceeds the number of samples. In this study, we analyzed ~28,600 genes and hundreds of lipids across 54 samples.

2. **Accommodates random effects**: ASCA implemented through the ALASCA package uses a linear mixed model framework, allowing inclusion of random effects (e.g., sample as a random effect) which is limited in traditional differential expression analyses like DESeq2.

3. **Captures coordinated responses**: Rather than examining each gene or lipid individually, ASCA identifies principal components that represent coordinated molecular responses across the entire dataset. This reveals underlying biological patterns that may be missed by univariate approaches.

## Benefits for Factorial Designs

4. **Partitions variance by experimental factors**: ASCA decomposes the total variance into contributions from main effects (temperature, parent) and interactions (temperature × parent), revealing which factors drive different patterns of molecular variation.

5. **Handles ordered categorical variables**: Temperature can be treated as an ordered categorical variable (27°C → 30°C → 33°C), allowing detection of linear and non-linear response patterns across the gradient.

6. **Robust statistical validation**: The ALASCA package implements bootstrap validation to assess the significance and reproducibility of identified components, ensuring patterns are not due to chance.

## Applicability to Multi-Omic Data

7. **Consistent framework across data types**: The same ASCA approach can be applied to different molecular data types (transcriptomics, metabolomics, lipidomics, proteomics) as long as the data structure is similar (samples × features matrix), enabling integrated multi-omic analyses.

8. **Appropriate for count-based and continuous data**: With proper normalization and transformation, ASCA handles both count-based data (RNAseq) and continuous measurements (lipid abundances).

# Analysis 1: RNAseq ASCA Analysis

**Script**: [`scripts/rna-seq_ASCA.Rmd`](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_ASCA.Rmd)

## Dataset and Experimental Design

The RNAseq dataset includes:
- **Samples**: 54 larval samples (6 biological replicates per treatment combination)
- **Genes**: ~28,600 genes after filtering (genes present at >10 counts in at least 10% of samples)
- **Design**: 3 temperatures × 3 parental phenotypes (full factorial)

## Analytical Approach

### Step 1: Data Preparation and Normalization

RNAseq count data requires appropriate normalization before ASCA analysis:

```r
# Filter genes using pOverA approach
filt <- filterfun(pOverA(0.1, 10))
gfilt <- genefilter(gcount, filt)
gkeep <- gcount[gfilt,]

# Apply variance stabilizing transformation (VST) using DESeq2
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                              colData = metadata_ordered,
                              design = ~temperature + parent + temperature:parent)

gvst <- vst(gdds, blind=TRUE)
normalized_counts <- assay(gvst)
```

**Why VST normalization?** Variance stabilizing transformation removes the mean-variance relationship in count data and stabilizes variance across the dynamic range, making the data suitable for PCA-based methods like ASCA. VST is preferred over log transformation for count data because it handles low counts better and provides more stable variance.

### Step 2: Format Data for ASCA

The normalized data is converted to long format required by the ALASCA package:

```r
long_data <- normalized_df %>%
  pivot_longer(cols = -c(sample, temperature, parent), 
               names_to = "variable", 
               values_to = "value")
```

### Step 3: ASCA Model Specification

The ASCA model is specified to test main effects and interactions:

```r
res_all <- ALASCA(
  long_data,
  value ~ temperature * parent + (1|sample),
  effects = c("temperature", "parent", "temperature:parent"),
  n_validation_runs = 1000,
  scale_function = "none",  # No additional scaling since data is VST-normalized
  validate = TRUE,
  stratification_column = "code",  # Ensures balanced resampling
  p_adjust_method = "fdr",
  reduce_dimensions = TRUE
)
```

**Key model features**:
- **Main effects**: Temperature and parent are modeled as fixed effects
- **Interaction term**: temperature:parent captures non-additive effects
- **Random effect**: (1|sample) accounts for sample-specific variation
- **No scaling**: Data is already variance-stabilized, so no additional scaling is needed
- **Stratified validation**: Ensures each bootstrap includes representatives from all treatment groups (critical with small sample sizes)
- **1000 validation runs**: Provides robust statistical assessment of component significance

### Step 4: Component Identification

ASCA performs PCA on the variance attributed to each effect separately. For each effect (temperature, parent, temperature:parent), the analysis identifies principal components and their variance explained. The first 10 components are examined for each effect.

### Step 5: Identifying Important Genes

For each principal component, genes are ranked by their loading scores (contribution to the component). "Important" genes are identified using a cumulative variance approach:

```r
# Rank genes by absolute loading
temp_pc1_load <- pc1_loadings %>%
  arrange(desc(abs_loading)) %>%
  mutate(
    var_contrib = loading^2,
    cum_var = cumsum(var_contrib) / sum(var_contrib)
  )

# Select genes explaining first 50% of variance
important_genes <- temp_pc1_load %>%
  filter(cum_var <= 0.50)
```

**Rationale for 50% threshold**: This stringent cutoff (tested against 60% and 80% thresholds) retains only the strongest drivers of each pattern while maintaining a manageable gene list size for functional enrichment. The number of genes identified (~3,000-4,000 per PC) is comparable to differential expression analyses but represents coordinated patterns rather than individual comparisons.

## Main Results

The ASCA analysis identified three major patterns of coordinated gene expression:

### Pattern 1 (Temperature PC1): Linear Temperature Response
- **Variance explained**: ~45%
- **Pattern**: Linear gradient across temperatures (27°C → 30°C → 33°C)
- **Genes identified**: 4,053 genes
- **Biological interpretation**: Core thermal stress response affecting nearly all larvae similarly regardless of parental origin

### Pattern 2 (Temperature:Parent PC1): Symbiont-Mediated Responses  
- **Variance explained**: ~20%
- **Pattern**: Separation between larvae with *Cladocopium*-only vs. mixed symbiont communities
- **Genes identified**: 3,247 genes
- **Biological interpretation**: Symbiont community identity drives baseline differences in gene expression, particularly pronounced at high stress (33°C)

### Pattern 3 (Temperature:Parent PC2): Parental History Effects
- **Variance explained**: ~15%
- **Pattern**: Unique expression profiles for each parental phenotype
- **Genes identified**: 3,225 genes  
- **Biological interpretation**: Parental history creates constitutive differences in offspring gene expression independent of temperature

**Minimal overlap between patterns**: An upset plot analysis showed that >90% of "important" genes were unique to each PC, confirming that these represent distinct biological patterns rather than redundant information.

## Computational Requirements

**Important**: This analysis requires substantial computational resources:
- **Memory**: Minimum 100 GB RAM
- **Runtime**: ~1.5 hours for full validation (1000 bootstrap runs)
- **Recommendation**: Run as a background job on HPC or server; use n=10-50 validation runs for exploratory analysis

# Analysis 2: Functional Enrichment of ASCA Results

**Script**: [`scripts/rna-seq_functional_enrichment_ASCA_topGO.Rmd`](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_functional_enrichment_ASCA_topGO.Rmd)

## Purpose

After identifying important genes for each biological pattern, functional enrichment analysis determines what biological processes and molecular functions these genes represent. This reveals the mechanisms underlying each pattern.

## Analytical Approach

### Step 1: Annotation Database Preparation

Functional annotations for *M. capitata* genes were obtained from the genome annotation (version 3):
- **Source**: EggNog functional annotation database
- **Genes annotated**: 24,072 genes with GO terms (from http://cyanophora.rutgers.edu/montipora/)
- **Coverage**: ~44% of predicted protein-coding genes in the genome
- **Detected in dataset**: 11,712 annotated genes present in our filtered RNAseq data

```r
# Load annotation file
Mcap.annot <- read.table("data/rna_seq/Montipora_capitata_HIv3.genes.EggNog_results.txt",
                         quote="", sep="\t", header=TRUE)

# Filter to only genes detected in our dataset
filtered_Mcap.annot <- Mcap.annot[Mcap.annot$gene %in% detected_genes, ]
```

### Step 2: TopGO Enrichment Analysis

TopGO (Topology-based Gene Ontology) analysis was performed for each PC gene list:

```r
# Prepare gene-to-GO mapping
geneID2GO <- filtered_Mcap.annot %>%
  select(gene, GOs) %>%
  filter(!is.na(GOs)) %>%
  separate_rows(GOs, sep = ";")

geneID2GO <- split(geneID2GO$GOs, geneID2GO$gene_id)

# Create binary vector (1 = important gene, 0 = background)
gene_vector <- as.factor(as.integer(background %in% important_gene_list))
names(gene_vector) <- background

# Build TopGO data object
GOdata <- new("topGOdata",
              description = "Temperature PC1",
              ontology = "BP",  # Biological Process
              allGenes = gene_vector,
              nodeSize = 5,  # Prune GO terms with <5 genes
              gene2GO = geneID2GO,
              annot = annFUN.gene2GO)

# Run weighted Fisher exact test
resultWeight <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

# Extract results
allRes <- GenTable(GOdata, weightFisher = resultWeight,
                  orderBy = "weightFisher",
                  topNodes = length(score(resultWeight)))
```

**Key parameters**:
- **nodeSize = 5**: Removes GO terms with fewer than 5 annotated genes to avoid spurious enrichments
- **algorithm = "weight01"**: Accounts for GO hierarchy topology, reducing redundancy
- **ontology = "BP"**: Focuses on Biological Process terms (also tested MF and CC separately)

### Step 3: Reducing GO Term Redundancy

GO terms are hierarchical and highly redundant. The rrvgo package reduces this redundancy:

```r
# Calculate GO term similarity
simMatrix <- calculateSimMatrix(significant_GO_terms,
                               orgdb="org.Ce.eg.db",
                               ont="BP",
                               method="Rel")

# Reduce to representative terms
reducedTerms <- reduceSimMatrix(simMatrix,
                               scores,
                               threshold=0.7,  # Similarity threshold
                               orgdb="org.Ce.eg.db")

# Group terms by parent
results_with_parents <- results %>%
  mutate(ParentTerm = reducedTerms$parentTerm[match(GO.ID, reducedTerms$go)])
```

**Effect of reduction**: The reduced list typically contains 100-150 specific GO terms grouped under 30-40 parent terms, making results more interpretable while preserving biological information.

### Step 4: Visualization

Results are visualized showing:
- **Gene Ratio**: Proportion of significant genes associated with each GO term
- **P-value**: Statistical significance (FDR-adjusted)
- **Parent terms**: Broader biological categories for interpretation

```r
# Calculate gene ratio
results$GeneRatio <- results$Significant / results$Annotated

# Plot parent terms
ggplot(parent_term_summary, 
       aes(x = GeneRatio, y = reorder(ParentTerm, GeneRatio))) +
  geom_point(aes(size = N_terms, color = -log10(min_pvalue))) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Mean Gene Ratio", y = "Parent GO Term",
       color = "-log10(P-value)", size = "N Terms")
```

## Main Results

Functional enrichment revealed distinct biological processes for each pattern:

### Temperature Response Pattern (PC1)
**Key processes** (4,053 genes enriched):
- Protein folding and unfolding (heat shock response)
- Oxidative stress response and DNA damage protection
- Cell cycle regulation and proliferation control
- Membrane trafficking and transport
- Metabolic reprogramming (shifts in energy pathways)
- Immune and symbiosis signaling
- Structural and developmental processes
- **Fatty acid biosynthesis** (relevant to lipid results below)

**Interpretation**: Classic thermal stress response showing cellular protection mechanisms, metabolic shifts, and stress signaling across all larvae.

### Symbiont-Mediated Pattern (PC2)
**Key processes** (3,247 genes enriched):
- Oxidative stress regulation and redox homeostasis
- Nitrogen transport and metabolism (ammonium handling)
- Metabolite transport and nutritional exchange
- Cell signaling and communication (host-symbiont crosstalk)
- Developmental and structural remodeling
- Cell cycle regulation

**Interpretation**: Constitutive differences in oxidative stress capacity, nitrogen management, and nutritional exchange between larvae with different symbiont communities. These functions relate to thermal tolerance differences observed in physiological assays.

### Parental History Pattern (PC3)
**Key processes** (3,225 genes enriched):
- Metabolic reprogramming (carbohydrate, fatty acid, amino acid metabolism)
- Stress responses and apoptosis regulation
- Development and cell structure
- Cell signaling and communication

**Interpretation**: Parental history creates baseline differences in metabolic programming and stress response capacity in offspring, suggesting transgenerational effects.

# Analysis 3: Lipidomics ASCA Analysis

**Script**: [`scripts/lipids_ASCA.Rmd`](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/lipids_ASCA.Rmd)

## Why ASCA is Appropriate for Lipidomic Data

ASCA is equally suitable for lipidomic data as for transcriptomic data because:

1. **Similar data structure**: Both datasets consist of a samples × features matrix with continuous measurements
2. **High dimensionality**: Lipidomic datasets typically include hundreds of lipid species, creating the same challenge of more features than samples
3. **Coordinated regulation**: Lipids within the same class or pathway often change together in response to stress, making them ideal for multivariate pattern detection
4. **Experimental design**: The same full factorial design (temperature × parent) applies to both data types
5. **Biological integration**: Lipid and gene expression changes should be coordinated (e.g., fatty acid biosynthesis genes and fatty acid lipids)

## Dataset and Experimental Design

The lipidomic dataset includes:
- **Samples**: Same 54 larval samples as RNAseq analysis
- **Lipid features**: Hundreds of distinct lipid species quantified
- **Design**: Same 3 temperatures × 3 parental phenotypes
- **Measurement**: Area-normalized peak intensities from MS-based lipidomics

## Analytical Approach

### Step 1: Data Preparation

Lipidomic data requires minimal preprocessing compared to RNAseq:

```r
# Load processed lipidomic data
lipids <- read_csv("output/lipids_metabolites/lipids/processed_lipids.csv")

# Convert to long format
long_data <- lipids %>%
  select(sample, temperature, parent, lipid, area_normalized) %>%
  rename(variable = lipid, value = area_normalized)

# Factor ordering (temperature as ordered categorical)
long_data$temperature <- factor(long_data$temperature,
                               levels=c("Ambient", "Moderate (+3)", "High (+6)"))
long_data$parent <- factor(long_data$parent,
                          levels=c("Wildtype", "Nonbleached", "Bleached"))
```

### Step 2: Data Transformation

Unlike VST-normalized RNAseq data, raw lipidomic abundances benefit from log transformation:

```r
# Check distribution
hist(long_data$value)  # Often right-skewed

# Apply log1p transformation (log(1 + x) to handle zeros)
long_data <- long_data %>%
  mutate(value = log1p(value))

hist(long_data$value)  # More normally distributed
```

**Why log1p transformation?** Lipidomic peak areas span several orders of magnitude with right-skewed distributions. Log transformation stabilizes variance and normalizes distributions, making the data more suitable for PCA. The log1p variant (log(1 + x)) handles zero values without requiring imputation.

### Step 3: ASCA Model Specification

The ASCA model for lipids is nearly identical to the RNAseq model:

```r
set.seed(123)

res_all <- ALASCA(
  long_data,
  value ~ temperature * parent + (1|sample),
  effects = c("temperature", "parent", "temperature:parent"),
  n_validation_runs = 1000,
  scale_function = "sdall",  # Scale by overall SD
  validate = TRUE,
  stratification_column = "code",
  p_adjust_method = "fdr",
  reduce_dimensions = TRUE
)
```

**Differences from RNAseq model**:
- **scale_function = "sdall"**: Scales all variables by overall standard deviation (appropriate for log-transformed data; RNAseq used "none" because VST already stabilizes variance)
- **Faster runtime**: ~5 minutes for full validation (vs. 1.5 hours for RNAseq) due to fewer features

### Step 4: Identifying Important Lipids

The same cumulative variance approach is used, but with a 25% threshold:

```r
# Rank lipids by loading contribution
temp_pc1_load <- pc1_loadings %>%
  arrange(desc(abs_loading)) %>%
  mutate(
    var_contrib = loading^2,
    cum_var = cumsum(var_contrib) / sum(var_contrib)
  )

# Find 25% cumulative variance threshold
idx_25_var <- which(temp_pc1_load$cum_var >= 0.25)[1]

# Select important lipids
important_lipids <- temp_pc1_load %>%
  filter(rank <= idx_25_var)
```

**Why 25% threshold for lipids?** Lipidomic datasets have fewer features than RNAseq (~hundreds vs. ~thousands), so a lower threshold still provides a reasonable number of important features while maintaining stringency.

## Main Results

ASCA identified similar patterns in lipidomic data as in gene expression:

### Pattern 1 (Temperature PC1): Linear Temperature Response in Lipids
- **Pattern**: Linear increase/decrease across temperature gradient
- **Biological interpretation**: Membrane lipid remodeling in response to thermal stress, consistent with temperature PC1 gene expression patterns including fatty acid biosynthesis genes

### Pattern 2 (Temperature PC2): Non-Linear Temperature Response
- **Pattern**: Higher/lower at moderate temperature relative to ambient and high
- **Biological interpretation**: Adaptive lipid remodeling at moderate stress that differs from extreme responses

### Pattern 3 (Parent PC1): Symbiont-Associated Lipid Differences
- **Pattern**: Higher/lower in Bleached vs. Wildtype and Nonbleached phenotypes
- **Biological interpretation**: Symbiont community identity affects lipid composition, consistent with gene expression PC2 showing symbiont-mediated effects

### Pattern 4 (Parent PC2): Parent-Specific Lipid Profiles
- **Pattern**: Higher/lower in Wildtype vs. other phenotypes
- **Biological interpretation**: Parental phenotype-specific lipid signatures

### Pattern 5 (Interaction PC1 & PC2): Temperature Response Varies by Parent
- **Pattern**: Responsive to temperature in Bleached and Nonbleached but not Wildtype
- **Biological interpretation**: Parental phenotype modulates the lipid remodeling response to temperature stress

## Integration with Gene Expression Results

The lipidomic ASCA results complement and validate the gene expression findings:

1. **Temperature responses**: Both datasets show strong linear temperature effects, with genes involved in fatty acid biosynthesis enriched in temperature PC1, consistent with lipid remodeling patterns
2. **Symbiont effects**: Both datasets separate by symbiont community, supporting the importance of symbiont identity in mediating molecular phenotypes
3. **Parental effects**: Both datasets show constitutive differences by parental history, supporting transgenerational effects on molecular phenotypes

# Comparison: RNAseq vs. Lipidomics ASCA

| Feature | RNAseq ASCA | Lipidomics ASCA |
|---------|-------------|-----------------|
| **Features** | ~28,600 genes | ~hundreds of lipids |
| **Normalization** | VST (DESeq2) | Log1p transformation |
| **Scaling** | None (pre-normalized) | sdall (scale by overall SD) |
| **Runtime** | ~1.5 hours | ~5 minutes |
| **Memory** | 100 GB | Standard (16-32 GB) |
| **Important features threshold** | 50% cumulative variance | 25% cumulative variance |
| **Functional enrichment** | TopGO (BP, MF, CC) | Lipid class enrichment |
| **Main patterns** | Temperature, symbiont, parent | Temperature, symbiont, parent |

# Methodological Considerations

## Strengths of ASCA for This Study

1. **Unified framework**: Same analysis approach for different -omic data types enables direct comparison
2. **Factorial design**: Explicitly models main effects and interactions rather than pairwise comparisons
3. **Coordinated patterns**: Identifies biological modules rather than individual features
4. **Random effects**: Accounts for sample-level variation not captured by treatment factors
5. **Validation**: Bootstrap approach provides robust assessment of pattern significance
6. **Ordered variables**: Treats temperature as ordered gradient rather than discrete groups

## Limitations and Considerations

1. **Normalization matters**: Different data types require appropriate transformation (VST for counts, log for abundances)
2. **Threshold selection**: Choice of cumulative variance threshold for "important" features is somewhat arbitrary (addressed by testing multiple thresholds)
3. **Computational cost**: RNAseq ASCA requires substantial computing resources for full validation
4. **Annotation coverage**: Functional enrichment limited to annotated genes (~44% of detected genes)
5. **Interpretation**: PC loadings indicate contribution but not directionality of effect (addressed by examining effect plots and example genes/lipids)

## Comparison to Alternative Approaches

**vs. DESeq2 differential expression**:
- ASCA captures coordinated patterns; DESeq2 identifies individual genes
- ASCA includes random effects; DESeq2 limited to fixed effects
- ASCA shows factor contributions; DESeq2 shows pairwise contrasts
- Both approaches are complementary rather than alternatives

**vs. WGCNA (Weighted Gene Co-expression Network Analysis)**:
- ASCA partitions by experimental factors; WGCNA clusters by correlation
- ASCA provides statistical validation; WGCNA identifies modules
- ASCA requires experimental design; WGCNA works with any dataset
- Both identify coordinated responses but from different perspectives

**vs. traditional MANOVA**:
- ASCA handles p >> n; MANOVA requires n >> p
- ASCA includes random effects; MANOVA limited to fixed effects
- ASCA provides PCA visualization; MANOVA provides test statistics
- ASCA is essentially a modern extension of MANOVA for high-dimensional data

# Key Takeaways

1. **ASCA is well-suited for full factorial designs in -omics studies**: The method explicitly models experimental factors and their interactions, providing clear interpretation of how temperature and parental effects structure molecular variation.

2. **Appropriate normalization is critical**: RNAseq requires variance-stabilizing transformation; lipidomics benefits from log transformation. The choice affects downstream scaling decisions.

3. **Coordinated patterns reveal biological mechanisms**: Rather than examining thousands of features individually, ASCA identifies 3-5 major patterns that represent distinct biological processes (thermal stress, symbiont effects, parental effects).

4. **Multi-omic integration is facilitated**: Using the same ASCA framework for different molecular data types enables direct comparison of patterns across data types (e.g., fatty acid biosynthesis genes and fatty acid lipids both respond linearly to temperature).

5. **Functional enrichment provides mechanistic insight**: Linking important genes to biological processes (TopGO) reveals what cellular mechanisms drive each pattern, moving from pattern description to biological interpretation.

6. **Validation ensures robustness**: Bootstrap validation with stratified resampling (1000 runs) provides confidence that patterns are biologically meaningful rather than statistical artifacts.

# Data and Code Availability

All analysis scripts and data are available in the project GitHub repository:

- **RNAseq ASCA script**: [`scripts/rna-seq_ASCA.Rmd`](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_ASCA.Rmd)
- **Functional enrichment script**: [`scripts/rna-seq_functional_enrichment_ASCA_topGO.Rmd`](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_functional_enrichment_ASCA_topGO.Rmd)
- **Lipidomics ASCA script**: [`scripts/lipids_ASCA.Rmd`](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/lipids_ASCA.Rmd)
- **RNAseq results**: [`output/rna_seq/asca/`](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/output/rna_seq/asca)
- **Lipidomics results**: [`output/lipids_metabolites/lipids/asca/`](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/output/lipids_metabolites/lipids/asca)
- **Figures**: [`figures/rna_seq/asca/`](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/figures/rna_seq/asca) and [`figures/lipids/asca/`](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/figures/lipids/asca)
- **Full repository**: [github.com/AHuffmyer/larval_symbiont_TPC](https://github.com/AHuffmyer/larval_symbiont_TPC)

# References

- Smilde, A.K. et al. (2005) ANOVA-simultaneous component analysis (ASCA): a new tool for analyzing designed metabolomics data. *Bioinformatics* 21(13):3043-3048
- Jansen, J.J. et al. (2005) ASCA: analysis of multivariate data obtained from an experimental design. *Journal of Chemometrics* 19:469-481
- Jarmund, A.H. et al. (2020) ALASCA: An R package for longitudinal and cross-sectional analysis of multivariate data by ASCA-based methods. *Frontiers in Molecular Biosciences* 7:585662
- Love, M.I. et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550
- Alexa, A. and Rahnenfuhrer, J. (2023) topGO: Enrichment Analysis for Gene Ontology. R package
- Trigg, S.A. et al. (2020) Transcriptomic responses to elevated pCO2 reveal threshold and cascade effects in an Antarctic pteropod. *BMC Genomics* 21:796
