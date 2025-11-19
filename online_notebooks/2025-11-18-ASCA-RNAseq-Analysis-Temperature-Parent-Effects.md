---
layout: post
title: ASCA Analysis of Temperature and Parent Effects on Gene Expression
date: '2025-11-18'
categories: Mcap_Larval_Symbiont_TPC
tags: Mcapitata RNAseq ASCA GeneExpression Temperature ParentalEffects R
---

This post details ANOVA Simultaneous Component Analysis (ASCA) of RNAseq data to examine how temperature and parental phenotype affect gene expression patterns in *Montipora capitata* larvae.

# Overview

In this analysis, I used ASCA (ANOVA Simultaneous Component Analysis) to examine coordinated gene expression patterns in coral larvae across temperatures (27°C, 30°C, and 33°C) and between parental phenotypes (Wildtype, Bleached, and Nonbleached). These phenotypes also vary in symbiont community composition with Wildtype and Nonbleached hosting a mixture of *Cladocopium* and *Durusdinium* with Bleached only hosting *Cladocopium*. ASCA (ALASCA package in R) is a multivariate method that combines ANOVA with PCA to identify the main patterns of variation in gene expression data and understand which genes drive those patterns.

Unlike differential expression analysis that examines genes one at a time, ASCA captures coordinated expression patterns (principle components) across the entire dataset, indicating how groups of genes respond together to experimental factors like temperature and phenotype effects in this study.

The data and analysis scripts for this project can be found in the [GitHub repository](https://github.com/AHuffmyer/larval_symbiont_TPC).

# Data and Methods

## Dataset

The analysis uses RNAseq data from *M. capitata* larvae that were:

- Exposed to three temperatures: 27°C (ambient), 30°C (moderate), and 33°C (high stress)
- Derived from three parental phenotypes: Wildtype (C and D symbionts), Bleached (parents with bleaching history; C symbints), and Nonbleached (parents without bleaching history; C and D symbints)  

The full dataset includes ~28,600 genes (after filtering) across 54 samples.  

## Analysis Approach

I performed ASCA using the ALASCA package in R, which uses a linear mixed model framework to:

1. **Normalize the data**: Applied variance stabilizing transformation (VST) to RNAseq counts using DESeq2
2. **Model the effects**: Used the model `value ~ temperature * parent + (1|sample)` to examine:
   - Main effect of temperature
   - Main effect of parent
   - Temperature × parent interaction
   - Sample as a random effect
3. **Identify components**: Used PCA to identify principal components that capture the most variation in gene expression
4. **Validate components**: Used bootstrapping (10 validation runs for exploratory analysis) to assess component significance
5. **Extract key genes**: Identified genes with the highest loadings (contributions) to each component (top 20)

The complete analysis script can be found at [`scripts/rna-seq_ASCA.Rmd`](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_ASCA.Rmd).  

# Results

## Scree Plot: Variance Explained

The scree plot shows how much variance is explained by each principal component across the first 10 components. We analyzed the first 10 components. 

![Scree Plot](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/asca/scree_plot.png)

**Interpretation**: 

The scree plot reveals that the first three principal components capture the majority of meaningful variation in the gene expression data:

- **PC1** explains approximately 40% of the variance, representing the dominant pattern of gene expression variation across temperature and parent treatments
- **PC2** explains approximately 15% of variance, capturing a secondary pattern
- **PC3** explains approximately 10% of variance, representing a third pattern of coordinated gene expression

Components beyond PC3 show diminishing returns with each explaining less than 5% of variance. This suggests that the main biological signals related to temperature and parent effects are captured in the first three components, which we'll examine in detail below. I've generated figures for PCs 4-6, but they explain only a minor portion of variance and are not reliable. 



## PC1: Primary Temperature Response Pattern

The first principal component captures the largest source of variation in gene expression and describes changes in gene expression linearly across temperature.

### PC1 Effect Plot

![PC1 Effect](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/asca/PC1.png)

**Interpretation**:

PC1 reveals clear separation between the temperature × parent treatment combinations:

- **Temperature gradient**: The component scores show a strong progression across temperatures (27°C → 30°C → 33°C), indicating that temperature is a major driver of PC1 variation
- **Parent effects**: Within each temperature, there are distinct differences between parental phenotypes
  - Wildtype larvae show intermediate responses
  - Bleached parent offspring show different expression patterns, potentially reflecting altered thermal tolerance
  - Nonbleached parent offspring show distinct patterns that may indicate enhanced resilience

The error bars (bootstrap confidence intervals) show that these differences are robust and reproducible. The clear separation suggests that PC1 captures a coordinated transcriptional response to thermal stress that differs by parental origin.

### PC1 Top Contributing Genes

![PC1 Genes](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/asca/PC1_genes.png)

**Interpretation**:

The genes contributing most to PC1 show coordinated expression patterns across treatments:

- **Thermal response genes**: The top genes show clear changes across the temperature gradient, with expression levels generally increasing or decreasing systematically from 27°C to 33°C
- **Parental modulation**: Several genes show differential responses based on parental phenotype, with bleached and nonbleached offspring showing divergent expression patterns compared to wildtype
- **Biological processes**: These genes likely represent core cellular responses to thermal stress, such as:
  - Heat shock proteins and molecular chaperones
  - Metabolic enzymes responding to changing energy demands
  - Stress response pathways
  - Cellular homeostasis maintenance genes

The coordinated changes in these genes drive the overall PC1 pattern, indicating they work together as part of an integrated thermal stress response that is modified by parental history.

## PC2: Secondary Response Pattern

The second principal component captures additional variation not explained by PC1.

### PC2 Effect Plot

![PC2 Effect](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/asca/PC2.png)

**Interpretation**:

PC2 reveals a more complex pattern than PC1:

- **Non-linear temperature response**: Unlike PC1's linear gradient, PC2 shows more variation at intermediate temperature (30°C), suggesting this component captures genes that respond differently at the middle temperature point
- **Enhanced parent differences**: The separation between parental phenotypes is more pronounced in PC2, particularly at 30°C and 33°C
- **Stress threshold effects**: The pattern suggests PC2 may capture genes that show threshold responses to thermal stress, where the response changes qualitatively at certain temperature points rather than gradually

This component likely represents secondary adaptation mechanisms or compensatory responses that become more important at elevated temperatures.

### PC2 Top Contributing Genes

![PC2 Genes](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/asca/PC2_genes.png)

**Interpretation**:

The genes driving PC2 show distinct patterns from PC1:

- **Variable responses**: These genes show more variability across treatments, with some increasing at intermediate temperatures and others showing peak expression at the extremes
- **Parent-specific patterns**: Several genes show strong expression differences between parental phenotypes, particularly distinguishing bleached from nonbleached offspring
- **Potential functions**: These genes may include:
  - Secondary stress response pathways
  - Parent-specific epigenetic programming
  - Symbiont-related genes (if symbionts are present)
  - Developmental regulation genes responding to thermal stress

The diversity of patterns in PC2 suggests it captures multiple biological processes that operate independently from the main thermal response in PC1.

## PC3: Tertiary Response Pattern

The third principal component captures additional independent variation.

### PC3 Effect Plot

![PC3 Effect](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/asca/PC3.png)

**Interpretation**:

PC3 shows yet another distinct pattern:

- **Complex temperature interaction**: The component scores show non-monotonic changes across temperature, with different parental groups showing distinct patterns
- **Parental phenotype divergence**: PC3 reveals strong separation between parental phenotypes, particularly at 33°C where the groups show maximal divergence
- **Potential acclimation signals**: The pattern suggests PC3 may capture genes involved in:
  - Longer-term acclimation responses
  - Parent-specific stress tolerance mechanisms
  - Trade-offs between different cellular processes

The larger confidence intervals in PC3 indicate more variability, which is expected for a tertiary component capturing more subtle or condition-specific responses.

### PC3 Top Contributing Genes

![PC3 Genes](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/asca/PC3_genes.png)

**Interpretation**:

The genes contributing to PC3 show:

- **Divergent expression patterns**: These genes show highly variable responses, with some showing peak expression at different temperatures depending on parental phenotype
- **Parent-specific responses**: Clear differences between how wildtype, bleached, and nonbleached offspring regulate these genes under thermal stress
- **Potential biological roles**:
  - Fine-tuning of stress responses
  - Parent-specific metabolic strategies
  - Developmental plasticity genes
  - Genes involved in long-term thermal acclimation

These genes represent more specialized responses that complement the major patterns captured in PC1 and PC2.

# Key Findings and Biological Insights

## Temperature Effects on Gene Expression

The ASCA analysis reveals that temperature has a profound and coordinated effect on gene expression:

1. **Dose-dependent response**: PC1 shows a clear temperature gradient from 27°C to 33°C, indicating that many genes respond proportionally to thermal stress
2. **Threshold effects**: PC2 and PC3 show non-linear patterns, suggesting some genes have temperature thresholds where their regulation changes qualitatively
3. **Coordinated regulation**: The high variance explained by PC1 indicates that hundreds of genes change expression together, representing an integrated cellular response to thermal stress

## Parental Effects on Gene Expression

Parental phenotype significantly modulates the thermal response:

1. **Persistent parent effects**: All three components show separation between offspring from different parental phenotypes, indicating that parental history affects offspring gene expression
2. **Differential thermal sensitivity**: Bleached and nonbleached offspring show distinct gene expression patterns compared to wildtype, suggesting different thermal tolerance strategies
3. **Interaction with temperature**: The changing separation between parental groups across temperatures (especially clear in PC2 and PC3) indicates that parent effects become more pronounced under thermal stress

## Transgenerational Effects

The results suggest potential transgenerational plasticity:

1. **Bleached offspring**: May carry molecular signatures of parental stress exposure, potentially affecting their own stress responses
2. **Nonbleached offspring**: Show distinct expression patterns that could indicate enhanced thermal tolerance inherited from resilient parents
3. **Epigenetic programming**: The parent-specific patterns in PC2 and PC3 suggest possible epigenetic modifications that persist to the offspring generation

## Biological Implications

These findings have important implications for coral reef resilience:

1. **Thermal tolerance variation**: Different parental backgrounds produce offspring with distinct molecular responses to thermal stress
2. **Potential for adaptation**: The variation in thermal responses across parental phenotypes suggests there is genetic and/or epigenetic diversity that selection could act upon
3. **Breeding programs**: Understanding parent-specific thermal response genes (as identified in this ASCA analysis) could inform selective breeding efforts for thermally tolerant corals
4. **Predictive markers**: The genes identified in PC1-3 could serve as biomarkers for thermal stress and resilience in *M. capitata*

# Next Steps

To build on this analysis, I plan to:

1. **Gene annotation**: Identify the biological functions of the top contributing genes in each component using GO term enrichment
2. **Network analysis**: Examine co-expression networks to understand how these genes interact
3. **Validation**: Confirm expression patterns of key genes using qPCR
4. **Integration**: Compare ASCA results with differential expression analysis to identify both coordinated patterns and individual gene changes
5. **Mechanism investigation**: Use the identified genes to develop hypotheses about molecular mechanisms underlying parental effects on thermal tolerance

# Data Availability

All data and analysis scripts are available on GitHub:

- **Analysis script**: [`scripts/rna-seq_ASCA.Rmd`](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_ASCA.Rmd)
- **Documentation**: [`scripts/rna-seq_ASCA_README.md`](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_ASCA_README.md)
- **Figures**: [`figures/rna_seq/asca/`](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/figures/rna_seq/asca)
- **Gene lists**: [`output/rna_seq/asca/`](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/output/rna_seq/asca)
- **Full repository**: [github.com/AHuffmyer/larval_symbiont_TPC](https://github.com/AHuffmyer/larval_symbiont_TPC)

# References

- Smilde, A.K. et al. (2005) ANOVA-simultaneous component analysis (ASCA): a new tool for analyzing designed metabolomics data. *Bioinformatics* 21(13):3043-3048
- Jansen, J.J. et al. (2005) ASCA: analysis of multivariate data obtained from an experimental design. *Journal of Chemometrics* 19:469-481
- Jarmund, A.H. et al. (2020) ALASCA: An R package for longitudinal and cross-sectional analysis of multivariate data by ASCA-based methods. *Frontiers in Molecular Biosciences* 7:585662
