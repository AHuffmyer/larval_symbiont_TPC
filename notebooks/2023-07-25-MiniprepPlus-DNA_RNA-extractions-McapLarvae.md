---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2023-07-25'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Larvae C/D 
---

# Extractions for Mcap larvae C/D Hawaii 2023 experiment 

This protocol is based on Maggie’s coral egg/larvae extraction [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/) and Kevin’s coral larvae extraction [protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/DNA-RNA-Extractions-on-P.-astreoides-larvae-BEAD-BEATING/). This post details the info about the extraction steps for the test samples that were collected in HI 2023 for A. Huffmyer's experiment on symbiont identity on *Montipora capitata* coral larval thermal performance. The github for that project is linked [here](https://github.com/AHuffmyer/larval_symbiont_TPC). 

### Samples 

The samples run today were R5, R8, and R41. These samples had 50 larvae in them and were preserved with 700 uL of DNA/RNA Shield. This photo is post-bead beating: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/samples_20230725.JPG)

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
- I eluted in 80 uL of either Tris (DNA) or DNA/RNA free water (RNA) instead of 100 uL. 

### QC 

#### Qubit & Nanodrop results 

| Sample ID | Qubit DNA1 | Qubit DNA2 | Qubit DNA average | Nanodrop DNA concentration | Nanodrop DNA 260/280 | Nanodrop DNA 260/230 | Qubit RNA1 | Qubit RNA2 | Qubit RNA average | Nanodrop RNA concentration | Nanodrop RNA 260/280 | Nanodrop RNA 260/230 |
| --------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- |
| R5        | 12.6       | 12.4       | 12.5              | 23.9                       | 1.94                 | 2.84                 | 19.4       | 19         | 19.2              | 9                          | 1.91                 | 1.19                 |
| R8        | 11.3       | 11.1       | 11.2              | 20.4                       | 1.87                 | 3.19                 | 21.4       | 21.2       | 21.3              | 9.8                        | 1.95                 | 0.66                 |
| R41       | 9.64       | 9.5        | 9.57              | 19.6                       | 1.96                 | 3.33                 | 17.8       | 17.2       | 17.5              | 8.6                        | 1.99                 | 1.04                 |

RNA concentrations look good! The 260/230 is still lower than I would like it to be. I could add another wash step or increase the volume of wash buffer (i.e., instead of washing in 400 uL, wash in 700 uL). 

#### Gel 

 ![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/gel_20230725.JPG)
 
 Gel looks good too! DNA still has those 18S and 28S bands that look like RNA, not sure why that is...