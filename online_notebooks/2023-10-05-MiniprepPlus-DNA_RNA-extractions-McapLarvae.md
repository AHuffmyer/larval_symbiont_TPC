---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2023-10-05'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Larvae C/D 
---

# Extractions for Mcap larvae C/D Hawaii 2023 experiment 

This protocol is based on Maggie’s coral egg/larvae extraction [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/) and Kevin’s coral larvae extraction [protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/DNA-RNA-Extractions-on-P.-astreoides-larvae-BEAD-BEATING/). This post details the info about the extraction steps for the test samples that were collected in HI 2023 for A. Huffmyer's experiment on symbiont identity on *Montipora capitata* coral larval thermal performance. The github for that project is linked [here](https://github.com/AHuffmyer/larval_symbiont_TPC). 

### Samples 

The samples run today were R20, R43, and R77. These samples had 50 larvae in them and were preserved with 700 uL of DNA/RNA Shield. This photo is post-bead beating: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/samples_20231005.JPG)

Today, I'm going to elute in 65 uL instead of 80 to see if I can concentrate the samples a little more. 

### Materials 

- Zymo Quick-DNA/RNA Miniprep Plus Kit [HERE](https://files.zymoresearch.com/protocols/_d7003t_d7003_quick-dna-rna_miniprep_plus_kit.pdf) Protocol Booklet
- Tris-Ethylenediaminetetraacetic acid (EDTA) 1X buffer for DNA elution
- Heating block capable of heating to 70ºC
- Centrifuge and rotor capable of spinning at 15,000 rcf
- Plastics 
	- 3 1.5 mL microcentrifuge tubes per sample
	- 2 PCR tubes per sample
	- 2 Qubit tubes per sample 
	- 1 5 mL tube per sample 

### Protocol

Protocol was followed according to this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-07-21-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae.md). A few differences: 

- I added 300 uL of DNA/RNA shield to the samples after thawing them so that the volume in the tubes was 1 mL of shield. I aliquoted out 500 uL for extractions and saved the remaining fraction in case we need to re-extract. 
- In the RNA extraction, I added another wash step with 400 uL of wash buffer after the DNase incubation (i.e., after the incubation, I did 400 uL of prep buffer, 700 uL of wash buffer, 400 uL of wash buffer, and another 400 uL of wash buffer).
- I eluted in 65 uL of either Tris (DNA) or DNA/RNA free water (RNA) instead of 100 uL. 

### QC 

#### Qubit & Nanodrop results 

| Sample ID | Qubit DNA1 | Qubit DNA2 | Qubit DNA average | Nanodrop DNA concentration | Nanodrop DNA 260/280 | Nanodrop DNA 260/230 | Qubit RNA1 | Qubit RNA2 | Qubit RNA average | Nanodrop RNA concentration | Nanodrop RNA 260/280 | Nanodrop RNA 260/230 |
| --------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- |
| R20       | 14.6       | 14.5       | 14.55             | 26.2                       | 2.01                 | 2.65                 | 13         | 13         | 13                | 11.3                       | 2                    | 1.13                 |
| R77       | 12.3       | 12.2       | 12.25             | 21.7                       | 1.92                 | 5.26                 | 12         | 11.8       | 11.9              | 8.9                        | 2.21                 | 1.28                 |
| R43       | 12.1       | 12         | 12.05             | 22                         | 1.9                  | 2.87                 | NA         | NA         | NA                | 8.5                        | 2.15                 | 1.13                 |

#### Gel 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/gel_20231005.JPG)

The concentrations are not as high as I expected when diluting in 65 uL. R43 also was too low for RNA on the Qubit side. 
