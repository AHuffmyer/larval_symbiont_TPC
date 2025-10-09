---
layout: post
title: Comprehensive Analysis Summary - Larval Symbiont TPC Project
date: '2025-10-09'
categories: Mcapitata_Larval_Symbiont_TPC
tags: Mcapitata Symbiodiniaceae TPC Respirometry RNA-seq ITS2 Metabolomics Lipidomics
---

This post provides a comprehensive summary of the analyses conducted in the *Montipora capitata* larval symbiont thermal performance curve (TPC) project from the Hawaii 2023 field season. This project investigated the effects of symbiont community identity on coral larval thermal tolerance through an integrated multi-omics approach.

# Project Overview

This research examined metabolic rates and molecular responses of *Montipora capitata* coral larvae across different temperature regimes and symbiont communities. The study was conducted at the Hawaii Institute of Marine Biology (HIMB) in summer 2023 and integrates physiological, molecular, and biochemical data to understand how symbiont identity influences larval thermal tolerance.

## Experimental Groups

Three main larval groups were examined:

- **Bleached (C)**: Larvae from bleached parent colonies (primarily *Cladocopium* symbionts)
- **Nonbleached (C/D)**: Larvae from resistant parent colonies (*Cladocopium* and *Durusdinium* symbionts)
- **Wildtype (C/D)**: Larvae from wild population (*Cladocopium* and *Durusdinium* symbionts)

## Temperature Treatments

- **Ambient**: 27°C
- **Moderate**: 30°C (+3°C)
- **High**: 33°C (+6°C)
- **Extreme**: 36°C (for some experiments)

---

# Respirometry Analyses

Multiple respirometry experiments were conducted using Sensor Dish Reader (SDR) systems to measure oxygen flux as a proxy for photosynthesis and respiration rates.

## 1. Thermal Performance Curves (TPC)

### Methods

Respiration rates were measured across a temperature gradient (27°C, 30°C, 33°C, 36°C, 39°C, 40°C) in dark conditions. Larvae (n=5 per well) were exposed to each temperature for 30-45 minutes in a 24-well SDR plate. Oxygen consumption was measured every 15 seconds.

**Analysis approach:**
- Rates extracted using LoLinR package with percentile rank method (alpha = 0.4)
- Blank corrections applied and rates normalized to larval number
- TPC curves fitted using `rTPC` package with multiple thermal performance models
- Model selection based on AIC values
- Thermal performance parameters estimated: Topt (thermal optimum), Tmax (critical thermal maximum), activation energy

**Scripts:**
- Rate extraction: `scripts/tpc_sdr_extraction.Rmd`
- Analysis and modeling: `scripts/tpc_sdr_analysis.Rmd`
- [View scripts on GitHub](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/scripts)

### Results

![TPC curves across symbiont types](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/tpc_sdr/TPCs.png)

**Figure 1.** Thermal performance curves showing respiration rates across temperature treatments for different symbiont groups. Data points represent individual measurements with fitted TPC curves overlaid.

![TPC parameter estimates](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/tpc_sdr/TPC_estimates.png)

**Figure 2.** Estimated thermal performance parameters including thermal optimum (Topt) and critical thermal maximum (Tmax) for each symbiont group.

**Key findings:**
- Temperature significantly affected respiration rates (ANOVA: p < 0.001)
- Thermal performance curves showed differences among symbiont groups
- Nonbleached larvae showed distinct thermal tolerance profiles
- Critical thermal maxima ranged from 34-37°C depending on symbiont type

---

## 2. Photosynthesis-Respiration (P-R) Measurements

### Methods

Photosynthesis and respiration were measured across a light gradient at ambient temperature (27°C). Light levels included:
- Dark (0 PAR) - dark respiration (Rd)
- Low light (50 PAR)
- Medium light (100 PAR)
- High light (500 PAR)
- Dark (0 PAR) - light enhanced respiration (LER)

Each phase lasted 15-20 minutes. From these measurements, multiple metabolic metrics were calculated:
- P net (net photosynthesis)
- P gross (gross photosynthesis = P net + Rd)
- P:R ratio (photosynthesis to respiration ratio)
- Light enhanced respiration (LER)

**Analysis approach:**
- Slopes extracted from each light phase using LoLinR
- Blank corrections and normalization to larval number
- Statistical comparisons using ANOVA with light × symbiont interactions
- Post-hoc Tukey HSD tests for pairwise comparisons

**Scripts:**
- Rate extraction: `scripts/larval_PR_sdr_extraction.Rmd`
- Analysis: `scripts/larval_PR_sdr_analysis.Rmd`

### Results

![Photosynthesis and respiration across light levels](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/pr_sdr/pr_dots.png)

**Figure 3.** Oxygen flux across light levels showing respiration (negative values) and photosynthesis (positive values) for each symbiont group. Note the increase in photosynthesis with increasing light intensity.

![P:R ratios](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/pr_sdr/pr_ratio.png)

**Figure 4.** Photosynthesis to respiration (P:R) ratios across symbiont groups, showing metabolic balance under different light conditions.

**Key findings:**
- Light significantly increased oxygen production (photosynthesis) across all groups (ANOVA: p < 0.001)
- Light enhanced respiration (LER) was higher than dark respiration (Rd) indicating symbiont-mediated metabolic stimulation
- P:R ratios showed autotrophic potential (>1) at medium and high light levels
- Symbiont type influenced photosynthetic capacity and P:R ratios

---

## 3. Photosynthesis-Irradiance (P-I) Curves

### Methods

Photosynthesis-irradiance (P-I) curves were generated by measuring photosynthesis across a gradient of light intensities (0, 50, 100, 200, 300, 400, 500 PAR) at ambient temperature. 

**Analysis approach:**
- P-I curves fitted to each well using non-linear regression
- Parameters estimated: 
  - Am (maximum photosynthetic rate)
  - AQY (alpha, photosynthetic efficiency)
  - Rd (dark respiration)
  - Ik (saturation irradiance)
  - Ic (compensation irradiance)
- Models fitted using `nls.multstart` package
- Statistical comparisons of P-I parameters across symbiont types

**Scripts:**
- Rate extraction: `scripts/larval_pi_curve_sdr_extraction.Rmd`
- Analysis: `scripts/larval_pi_curve_sdr_analysis.Rmd`

### Results

![P-I curves](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/larval_pi_curves/nls_curves_all.png)

**Figure 5.** Photosynthesis-irradiance curves showing photosynthetic responses to increasing light intensity. Curves are fitted using non-linear regression with estimated parameters.

![P-I parameters](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/larval_pi_curves/pi_means.png)

**Figure 6.** Mean photosynthetic rates across light levels for each symbiont group, showing standard error bars.

**Key findings:**
- Photosynthesis saturated at approximately 100-200 PAR
- Maximum photosynthetic capacity (Am) differed by symbiont type
- Photosynthetic efficiency (AQY) was similar across groups
- Saturation irradiance (Ik) ranged from 50-100 PAR

---

## 4. Phenoplate Photophysiology Measurements

### Methods

Photophysiology was assessed using a Phenoplate fluorometer measuring:
- Quantum yield (QY) - photosynthetic efficiency
- Non-photochemical quenching (NPQ) - photoprotection
- Rapid light curves (RLC) - photosynthetic performance across light intensities

Measurements were conducted across temperature treatments (26°C, 30°C, 33°C) and symbiont groups.

**Analysis approach:**
- QY and NPQ calculated from fluorescence measurements
- RLC parameters fitted: rETRmax (maximum electron transport rate), alpha (photosynthetic efficiency), Ik (saturation irradiance)
- Statistical analyses using ANOVA with temperature × symbiont interactions
- Phase-specific analyses for light adaptation, dark adaptation, and stress responses

**Scripts:**
- RLC analysis: `scripts/phenoplate_RLC_analysis.Rmd`
- Light-dark phase analysis: `scripts/phenoplate_light_dark_phases_analysis.Rmd`

### Results

![RLC parameters](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/phenoplate/rlc_rETRmax.jpeg)

**Figure 7.** Maximum relative electron transport rate (rETRmax) across symbiont groups, showing photosynthetic capacity.

![NPQ responses](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/phenoplate/rlc_npq_symbiont.jpeg)

**Figure 8.** Non-photochemical quenching (NPQ) responses showing photoprotective mechanisms across symbiont types.

**Key findings:**
- rETRmax decreased with increasing temperature stress
- NPQ increased under high light and high temperature conditions
- Symbiont identity affected photophysiological responses
- Nonbleached larvae showed enhanced photoprotection

---

## 5. F1/F2 Selective Breeding Study

### Methods

A parallel study examined larvae from selectively bred coral lineages in collaboration with the Coral Resilience Lab at HIMB. Groups included:
- F1 Wildtype (WT) - from wild population
- F1 Non-bleached (NB) - from bleaching-resistant parents
- F2 WT - second generation wildtype
- F2 NB - second generation non-bleached

Respirometry was conducted at 27°C, 33°C, and 36°C across the same light gradient as the main P-R experiments.

**Analysis approach:**
- Same methodology as main P-R experiments
- ANOVA with temperature × group × light interactions
- Comparisons between F1 and F2 generations and between breeding lines

**Scripts:**
- Rate extraction: `scripts/F1F2_sdr_extraction.Rmd`
- Analysis: `scripts/F1F2_sdr_analysis.Rmd`

More details can be found in the [dedicated F1F2 project notebook post](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/online_notebooks/2023-07-31-F1F2-project-data-analysis-Hawaii-2023.md).

---

# Molecular Analyses

## 1. ITS2 Symbiont Community Analysis

### Methods

Symbiont community composition was characterized using ITS2 amplicon sequencing of the ribosomal internal transcribed spacer 2 region. 

**Sequencing and bioinformatics:**
- DNA extracted from larval samples (n=~60)
- ITS2 region amplified using universal primers
- Illumina sequencing (2x250bp paired-end)
- Quality control with FastQC/MultiQC
- Sequences analyzed using SymPortal framework
- ITS2 type profiles and DIVs (distinct ITS2 variants) identified

**Statistical analysis:**
- Alpha diversity (Shannon index, richness)
- Beta diversity using NMDS ordination (Bray-Curtis dissimilarity)
- PERMANOVA to test for group differences
- Heatmaps showing relative abundance of ITS2 profiles
- Community composition visualized by profile and DIV level

**Scripts:**
- Analysis: `scripts/its2.Rmd`
- [SymPortal outputs](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/data/its2/symportal)

### Results

![NMDS ordination](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/its2/nmds_its2_div.jpeg)

**Figure 9.** Non-metric multidimensional scaling (NMDS) ordination of symbiont communities based on ITS2 DIV profiles. Points colored by parent colony type show clear separation of symbiont communities.

![ITS2 profile composition](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/its2/profile_plot_means.png)

**Figure 10.** Relative abundance of dominant ITS2 type profiles across larval groups. Bleached larvae (C) are dominated by *Cladocopium* (C-type) while nonbleached and wildtype larvae contain both *Cladocopium* and *Durusdinium* (D-type) symbionts.

![Alpha diversity](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/its2/alpha_rare.png)

**Figure 11.** Alpha diversity rarefaction curves showing symbiont diversity within samples.

**Key findings:**
- Bleached larvae harbored primarily *Cladocopium* (genus C)
- Nonbleached and wildtype larvae contained mixed *Cladocopium* and *Durusdinium* (genus D) communities
- Symbiont communities were significantly different among parent colony types (PERMANOVA: p < 0.001)
- *Durusdinium* presence correlated with bleaching resistance history
- Multiple ITS2 type profiles detected within each symbiont genus

---

## 2. RNA-Sequencing and Differential Gene Expression

### Methods

Transcriptomic responses to thermal stress were assessed using RNA-sequencing.

**Library preparation and sequencing:**
- Total RNA extracted from larvae (n=48 samples)
- Quality assessed via TapeStation and Qubit
- RNA-seq libraries prepared (polyA selection)
- Illumina NovaSeq sequencing (150bp paired-end)
- ~30-40 million reads per sample

**Bioinformatics pipeline:**
- Quality control: FastQC, MultiQC
- Adapter trimming: fastp
- Alignment to *M. capitata* genome (HIv3) using HISAT2
- Gene quantification with featureCounts
- Detailed pipeline: `scripts/bioinformatics/2024-02-23-RNAseq-QC-for-Mcapitata-Larval-Thermal-Tolerance-Project.md`

**Differential expression analysis:**
- DESeq2 package in R
- Design: ~parent + temperature + parent:temperature
- Multiple contrasts tested:
  - Parent effects (Bleached vs Nonbleached vs Wildtype)
  - Temperature effects (27°C vs 30°C vs 33°C)
  - Interaction effects (parent × temperature)
- Significance threshold: padj < 0.05, |log2FC| > 1
- PCA for sample clustering and quality control
- Heatmaps of differentially expressed genes (DEGs)
- Volcano plots showing magnitude and significance of changes

**Scripts:**
- Differential expression: `scripts/rna-seq_DEG.Rmd`
- Functional enrichment: `scripts/rna-seq_functional_enrichment_goseq.Rmd`
- GO analysis: `scripts/rna-seq_functional_enrichment_topGO.Rmd`

### Results

![PCA plot](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/pca.jpeg)

**Figure 12.** Principal component analysis (PCA) of gene expression showing sample clustering by temperature and parent colony type.

![Volcano plot - parent effects](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/volcano_parent.png)

**Figure 13.** Volcano plot showing differentially expressed genes between parent colony types. Points above the horizontal line meet significance threshold (padj < 0.05) and points outside vertical lines meet fold-change threshold (|log2FC| > 1).

![Volcano plot - temperature effects](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/volcano_temperature.png)

**Figure 14.** Volcano plot showing differentially expressed genes across temperature treatments.

![PCA of DEGs - parent](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/pca_DEG_parent.png)

**Figure 15.** PCA using only differentially expressed genes for parent effects, showing clear separation of groups.

![PCA of DEGs - temperature](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/rna_seq/pca_DEG_temperature.png)

**Figure 16.** PCA using only differentially expressed genes for temperature effects.

**Key findings:**
- **Parent effects**: 1,500+ genes differentially expressed between parent types
  - Genes related to stress response, immune function, and metabolism enriched
  - Nonbleached larvae showed elevated expression of heat shock proteins (HSPs)
  - Cell cycle and development genes differed among groups
  
- **Temperature effects**: 2,000+ genes responsive to thermal stress
  - Heat shock response activated at 30°C and 33°C
  - Metabolic genes downregulated at high temperature
  - Oxidative stress response genes upregulated
  - Apoptosis and cell death pathways activated at 33°C

- **Interaction effects**: Parent-specific temperature responses detected
  - Nonbleached larvae showed attenuated stress response at moderate temperature
  - Different parent types employed distinct molecular strategies for thermal tolerance

**Functional enrichment:**
- GO terms enriched in temperature-responsive genes included:
  - Protein folding and unfolded protein response
  - Response to heat and oxidative stress
  - Metabolic processes
  - Apoptosis and programmed cell death
  - Immune response

---

## 3. Metabolomics

### Methods

Polar metabolites were profiled using gas chromatography-mass spectrometry (GC-MS).

**Sample preparation:**
- Metabolites extracted from larvae using methanol-chloroform-water extraction
- Polar phase collected and derivatized for GC-MS
- Internal standards added for normalization

**Analysis:**
- GC-MS analysis (compound identification and quantification)
- Peak integration and normalization to protein content
- Data processed using MS-Dial

**Statistical analysis:**
- Multivariate analysis: PCA, PLS-DA (Partial Least Squares Discriminant Analysis)
- Univariate analysis: ANOVA for each metabolite
- VIP scores (Variable Importance in Projection) to identify discriminating metabolites
- Pathway analysis linking metabolites to biological processes

**Scripts:**
- Analysis: `scripts/metabolomics.Rmd`

### Results

**Key findings:**
- Metabolite profiles differed significantly by parent type and temperature
- Amino acid levels changed under thermal stress
- Energy metabolism (TCA cycle intermediates) affected by temperature
- Osmolyte concentrations varied among groups
- Oxidative stress markers elevated at high temperature

---

## 4. Lipidomics

### Methods

Lipid composition was analyzed to understand membrane remodeling and energy storage responses.

**Sample preparation:**
- Lipids extracted from non-polar phase of metabolite extractions
- Analyzed by liquid chromatography-mass spectrometry (LC-MS)

**Analysis:**
- Lipid identification and quantification
- Normalized to protein content
- Lipid classes and individual lipid species characterized

**Statistical analysis:**
- PCA and PLS-DA for multivariate patterns
- VIP scores to identify important lipids
- Class-level analysis (phospholipids, glycerolipids, etc.)
- Specific fatty acid analysis (EPA, DHA, omega-3/omega-6 ratios)
- ANOVA for individual lipids and lipid classes

**Scripts:**
- Analysis: `scripts/lipidomics_area.Rmd`
- Protein normalization: `scripts/protein_omic_extractions.Rmd`

### Results

![Lipid PCA](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/lipids/pca_all_groups.jpeg)

**Figure 17.** Principal component analysis of lipid profiles showing separation by parent colony type and temperature treatment.

![Total lipids](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/lipids/total_lipids.jpeg)

**Figure 18.** Total lipid content across experimental groups, normalized to protein content.

![VIP scores - parent](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/lipids/VIP_scores_parent.png)

**Figure 19.** Variable Importance in Projection (VIP) scores identifying lipids that discriminate between parent types.

![VIP scores - temperature](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/lipids/VIP_scores_temperature.png)

**Figure 20.** VIP scores for temperature effects on lipid profiles.

![EPA/DHA ratio](https://raw.githubusercontent.com/AHuffmyer/larval_symbiont_TPC/main/figures/lipids/epa_dha_ratio.png)

**Figure 21.** Ratio of eicosapentaenoic acid (EPA) to docosahexaenoic acid (DHA), important omega-3 fatty acids, across groups.

**Key findings:**
- Total lipid content varied by parent type and temperature
- Membrane lipid remodeling occurred under thermal stress
- Phospholipid composition changed at high temperature
- Polyunsaturated fatty acids (PUFAs) decreased at elevated temperature
- EPA and DHA levels affected by symbiont type
- Nonbleached larvae maintained higher PUFA levels under stress
- Lipid saturation increased with temperature (membrane rigidification)

---

# Environmental and Morphological Analyses

## 1. Environmental Monitoring

### Methods

Environmental conditions were monitored throughout experiments:
- Temperature (loggers and daily measurements)
- Light intensity (PAR sensors)
- pH (total and NBS scales)
- Salinity
- Flow rates

**Analysis approach:**
- Time series analysis of logger data
- Statistical comparisons of daily measurements among treatment tanks
- ANOVA to verify treatment consistency
- Quality control checks

**Scripts:**
- Logger analysis: `scripts/logger_analysis.Rmd`
- Daily measurements: `scripts/daily_measurement_analysis.Rmd`

### Results

**Key findings:**
- Temperature treatments maintained target levels ± 0.5°C
- Light levels consistent across treatment tanks (p > 0.05)
- pH, salinity, and flow rates not significantly different among tanks
- Environmental conditions stable throughout experiments

---

## 2. Cell Density and Size Measurements

### Methods

Symbiont cell density and size were quantified to understand symbiont physiology.

**Measurements:**
- Cell counts using hemocytometer
- Cell size measured from microscopy images
- Normalized to larval number

**Analysis approach:**
- ANOVA comparing cell density and size across groups
- Correlations with metabolic rates
- Visualization of cell density distributions

**Scripts:**
- Analysis: `scripts/cells_size.Rmd`

### Results

**Key findings:**
- Cell density varied among symbiont types
- Cell size consistent within symbiont genera
- Cell density correlated with photosynthetic capacity
- Temperature affected cell density in some groups

---

# Summary and Integration

This comprehensive multi-omics study reveals:

## Physiological Responses
- Respiration rates increased with temperature up to thermal optima (~30-33°C)
- Critical thermal maxima differed by symbiont type
- Photosynthetic capacity influenced by symbiont identity
- Nonbleached larvae showed enhanced thermal tolerance in some metrics

## Molecular Responses
- Thousands of genes responsive to thermal stress
- Parent-specific gene expression patterns
- Heat shock response and metabolic reprogramming under stress
- Nonbleached larvae showed distinct molecular signatures

## Biochemical Responses
- Metabolite profiles reflected stress states
- Lipid remodeling occurred under thermal stress
- Membrane composition changes related to temperature adaptation
- Energy metabolism affected by temperature and symbiont type

## Integration
- Symbiont community composition (*Cladocopium* vs *Cladocopium*/*Durusdinium* mixes) influenced physiological, molecular, and biochemical responses
- Presence of *Durusdinium* associated with enhanced thermal tolerance in some contexts
- Multi-level responses (gene expression → metabolism → physiology) showed coordinated stress responses
- Parent bleaching history influenced offspring responses, suggesting potential for adaptation

---

# Data and Code Availability

All data and analysis scripts are available on GitHub:
- **Repository**: [https://github.com/AHuffmyer/larval_symbiont_TPC](https://github.com/AHuffmyer/larval_symbiont_TPC)
- **Scripts**: [https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/scripts](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/scripts)
- **Figures**: [https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/figures](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/figures)

RNA-seq and ITS2 sequence data have been deposited in NCBI SRA.

---

# Acknowledgments

This work was conducted in collaboration with the Coral Resilience Lab at the Hawaii Institute of Marine Biology. Thanks to the Putnam Lab at the University of Rhode Island for support and resources.

**Project contributors:**
- Ariana S. Huffmyer (analysis, data collection)
- Putnam Lab, University of Rhode Island
- Coral Resilience Lab, Hawaii Institute of Marine Biology

For more information, visit:
- [Putnam Lab GitHub](https://github.com/Putnam-Lab)
- [Coral Resilience Lab](https://www.coralresiliencelab.com/)

---

*This notebook post summarizes analyses conducted on data collected during the summer 2023 field season at HIMB. For detailed methods and results, please refer to the individual analysis scripts and the project README in the GitHub repository.*
