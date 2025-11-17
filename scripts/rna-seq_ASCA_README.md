# RNAseq ASCA Analysis

## Overview

This script (`rna-seq_ASCA.Rmd`) performs ANOVA Simultaneous Component Analysis (ASCA) on RNAseq data to examine how gene expression changes across temperatures (27, 30, and 33°C) for each parental phenotype (wildtype, bleached, and nonbleached).

## Purpose

ASCA is a multivariate statistical method that combines ANOVA with Principal Component Analysis (PCA) to:

- Decompose variation in gene expression according to experimental factors (temperature, parent phenotype)
- Identify the main patterns of gene expression changes across temperatures
- Determine which genes contribute most to thermal responses
- Reveal phenotype-specific thermal response signatures
- Identify thermally responsive genes unique to each parental phenotype

Unlike traditional differential expression analysis that tests genes one-at-a-time, ASCA captures coordinated patterns across the entire transcriptome.

## Method

The script uses **ALASCA (ANOVA Linear ASCA)** to:
1. Apply variance stabilizing transformation (VST) to normalize gene counts
2. Separate the analysis by parental phenotype (wildtype, bleached, nonbleached)
3. Model temperature effects on gene expression for each phenotype
4. Use bootstrapping (100 validation runs) to estimate component significance
5. Identify principal components explaining temperature-related variation
6. Extract genes with highest loadings (contributions) to each component
7. Compare thermally responsive genes across phenotypes to find unique signatures

## Requirements

### Input Files
- `data/rna_seq/Mcapitata2023_gene_count_matrix_STAR.csv` - Gene count matrix from STAR alignment
- `data/rna_seq/sample_rnaseq_metadata.csv` - Sample metadata with temperature and parent information

### R Packages
Install required packages:
```r
# BiocManager packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("genefilter", "DESeq2"))

# CRAN packages
install.packages(c("tidyverse", "RColorBrewer"))

# ALASCA from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("andjar/ALASCA", ref = "main")
```

## Analysis Workflow

### 1. Data Preparation
- Loads gene count matrix and metadata
- Removes genes with zero counts across all samples
- Applies pOverA filtering (genes must have >10 counts in 10% of samples)
- Reduces ~54,000 genes to ~28,600 genes with sufficient expression

### 2. Normalization
- Creates DESeq2 dataset with temperature + parent + interaction design
- Calculates size factors to normalize for sequencing depth
- Applies Variance Stabilizing Transformation (VST)
  - Reduces mean-variance relationship
  - Makes data approximately homoskedastic
  - Suitable for multivariate analysis

### 3. ASCA Analysis by Phenotype

The analysis is performed separately for each parental phenotype:

#### Wildtype
- Filters data to wildtype parents only
- Models: `value ~ temperature + (1|sample)`
- Temperature as fixed effect, sample as random effect
- Identifies components explaining thermal response variation

#### Bleached
- Filters data to bleached parents only
- Same model structure as wildtype
- Reveals thermal response in corals from bleached parents

#### Nonbleached
- Filters data to nonbleached parents only
- Same model structure
- Shows thermal response in corals from healthy parents

### 4. Component Analysis

For each phenotype:
- **Scree plot**: Shows variance explained by each principal component
- **Effect plots**: Display how component scores change across temperatures
- **Gene plots**: Show expression patterns of top contributing genes
- **Loading extraction**: Identifies top 20 genes for each component

### 5. Comparative Analysis

- Combines gene lists from all components per phenotype
- Identifies unique thermally responsive genes for each phenotype
- Creates summary showing phenotype-specific thermal signatures

## Outputs

### Figures
All saved to `output/rna_seq/asca/[phenotype]/`:

**Scree Plots:**
- `scree_wildtype.pdf` - Variance explained by components (wildtype)
- `scree_bleached.pdf` - Variance explained by components (bleached)
- `scree_nonbleached.pdf` - Variance explained by components (nonbleached)

**Effect Plots (PC1-3 for each phenotype):**
- `PC1_wildtype.pdf` - Component 1 scores across temperatures
- `PC2_wildtype.pdf` - Component 2 scores across temperatures
- `PC3_wildtype.pdf` - Component 3 scores across temperatures
- (Similar files for bleached and nonbleached)

**Gene Plots (PC1-3 for each phenotype):**
- `PC1_genes_wildtype.pdf` - Top 10 genes for component 1
- `PC2_genes_wildtype.pdf` - Top 10 genes for component 2
- `PC3_genes_wildtype.pdf` - Top 10 genes for component 3
- (Similar files for bleached and nonbleached)

### Data Files
All saved to `output/rna_seq/asca/`:

**Gene Lists per Phenotype:**
- `wildtype/wildtype_thermally_responsive_genes.csv` - Top 20 genes per component
- `bleached/bleached_thermally_responsive_genes.csv` - Top 20 genes per component
- `nonbleached/nonbleached_thermally_responsive_genes.csv` - Top 20 genes per component

**Comparative Analysis:**
- `unique_thermally_responsive_genes_by_phenotype.csv` - Genes unique to each phenotype

## Interpretation

### What the Results Tell You

1. **Scree Plots**: 
   - Components with >1-5% variance are typically meaningful
   - First few components usually capture most thermal variation
   - Differences in component structure between phenotypes indicate distinct response strategies

2. **Effect Plots**:
   - Show how component scores change across 27, 30, and 33°C
   - Linear patterns suggest gradual temperature response
   - Non-linear patterns indicate threshold effects or acclimation limits
   - Ribbons show bootstrap confidence intervals

3. **Gene Plots**:
   - Display normalized expression of top contributing genes
   - Genes with high loadings drive the temperature response
   - Similar patterns across genes indicate coordinated regulation
   - Different patterns between phenotypes reveal distinct mechanisms

4. **Unique Genes**:
   - Genes appearing only in one phenotype's top lists
   - Represent phenotype-specific thermal adaptation strategies
   - Potential genetic markers for thermal tolerance
   - Candidates for follow-up functional studies

5. **Component Loadings**:
   - Positive loadings: genes increasing with temperature
   - Negative loadings: genes decreasing with temperature
   - Magnitude indicates strength of contribution

### Biological Insights

- **Wildtype**: Baseline thermal response patterns
- **Bleached**: May show stress response signatures or altered resilience
- **Nonbleached**: May show enhanced thermal tolerance mechanisms
- **Unique genes**: Identify parental effects on offspring thermal performance

## Running the Script

1. Ensure input files are in `data/rna_seq/` directory
2. Install required R packages (see above)
3. Open `rna-seq_ASCA.Rmd` in RStudio
4. Run chunks sequentially or knit the entire document
5. Check `output/rna_seq/asca/` for results

## Computational Considerations

- **Runtime**: Analysis may take 60-90 minutes depending on your system
  - Most time spent on ASCA bootstrapping (100 validation runs per phenotype)
  - Three separate ASCA models (one per phenotype)
- **Memory**: Requires ~8-16 GB RAM
  - Full gene expression matrix is large (~28,600 genes × 54 samples)
- **Parallelization**: ALASCA can use multiple cores if available

## Customization

You can modify the script to:

1. **Change validation runs**:
   - Increase `n_validation_runs` (currently 100) for more robust estimates
   - Decrease for faster computation during testing

2. **Adjust gene filtering**:
   - Modify pOverA threshold (currently 0.1, 10)
   - More stringent filtering reduces genes but improves signal-to-noise

3. **Extract more genes**:
   - Change `n_limit` parameter (currently 20 for lists, 10 for plots)
   - Useful for finding additional candidate genes

4. **Analyze additional components**:
   - Add PC4, PC5 plots for components with >1% variance
   - May reveal subtle thermal response patterns

5. **Combined phenotype analysis**:
   - Could add analysis with both temperature and parent as main effects
   - Would show overall patterns but lose phenotype-specific resolution

## Differences from DEG Analysis

**ASCA Analysis (this script)**:
- Multivariate approach considering all genes together
- Identifies coordinated expression patterns
- Captures temperature gradients (27, 30, 33°C)
- Focuses on continuous variation
- Reveals phenotype-specific thermal signatures

**DEG Analysis** (`rna-seq_DEG.Rmd`):
- Univariate approach testing each gene separately
- Identifies significantly differentially expressed genes
- Typically compares discrete conditions (e.g., 27 vs 33°C)
- Focuses on statistical significance
- Provides fold changes and p-values

**Both approaches are complementary**: ASCA finds patterns, DEG finds specific changes.

## References

- Smilde, A.K. et al. (2005) ANOVA-simultaneous component analysis (ASCA): a new tool for analyzing designed metabolomics data. *Bioinformatics* 21(13):3043-3048
- Jansen, J.J. et al. (2005) ASCA: analysis of multivariate data obtained from an experimental design. *Journal of Chemometrics* 19:469-481
- Jarmund, A.H. et al. (2020) ALASCA: An R package for longitudinal and cross-sectional analysis of multivariate data by ASCA-based methods. *Frontiers in Molecular Biosciences* 7:585662

## Contact

For questions about this analysis, contact the Huffmyer Lab or open an issue in the repository.
