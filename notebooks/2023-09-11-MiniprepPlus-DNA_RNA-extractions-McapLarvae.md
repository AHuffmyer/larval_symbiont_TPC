---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2023-09-11'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Larvae C/D 
---

# Extractions for Mcap larvae C/D Hawaii 2023 experiment 

This protocol is based on Maggie’s coral egg/larvae extraction [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/) and Kevin’s coral larvae extraction [protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/DNA-RNA-Extractions-on-P.-astreoides-larvae-BEAD-BEATING/). This post details the info about the extraction steps for the test samples that were collected in HI 2023 for A. Huffmyer's experiment on symbiont identity on *Montipora capitata* coral larval thermal performance. The github for that project is linked [here](https://github.com/AHuffmyer/larval_symbiont_TPC). 

### Samples 

The samples run today were R1, R10, R17, R23, R29, R39, R44, R53, R59, R69, R79, R94, R99, R107. These samples had 50 larvae in them and were preserved with 700 uL of DNA/RNA Shield. This photo is post-bead beating: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/samples_20230911.JPG)


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
| R01       | 14.3       | 14.6       | 14.45             | 22.2                       | 1.99                 | 1.92                 | 11.8       | 11.6       | 11.7              | 10.3                       | 1.98                 | 1.31                 |
| R10       | 6.6        | 6.78       | 6.69              | 12.2                       | 1.84                 | 1.17                 | NA         | NA         | NA                | 7.7                        | 2.02                 | 1.11                 |
| R17       | 8.64       | 8.86       | 8.75              | 17.9                       | 2.01                 | 1.72                 | NA         | NA         | NA                | 10.5                       | 2.05                 | 0.1                  |
| R23       | 9.16       | 9.34       | 9.25              | 23.5                       | 1.95                 | 1.87                 | NA         | NA         | NA                | 9.2                        | 2.09                 | 1.21                 |
| R29       | 15.1       | 15.4       | 15.25             | 21.3                       | 1.97                 | 1.32                 | NA         | NA         | NA                | 8.8                        | 2.18                 | 1.14                 |
| R39       | 16.6       | 17         | 16.8              | 26                         | 1.96                 | 1.88                 | 13.8       | 13.8       | 13.8              | 13.5                       | 2.05                 | 1.27                 |
| R44       | 10.3       | 10.5       | 10.4              | 15.6                       | 1.92                 | 1.95                 | NA         | NA         | NA                | 9.6                        | 2.08                 | 1.18                 |
| R53       | 7.18       | 7.28       | 7.23              | 14.6                       | 2.01                 | 1.89                 | NA         | NA         | NA                | 7.3                        | 2.08                 | 1.06                 |
| R59       | 11.2       | 11.2       | 11.2              | 28.2                       | 1.93                 | 0.13                 | NA         | NA         | NA                | 7.1                        | 2.13                 | 1.2                  |
| R69       | 13.9       | 14.3       | 14.1              | 14.9                       | 1.9                  | 1.91                 | NA         | NA         | NA                | 7.7                        | 2                    | 1.22                 |
| R79       | 35.4       | 36         | 35.7              | 47.1                       | 1.93                 | 2.15                 | 16.8       | 16.4       | 16.6              | 16.9                       | 2.05                 | 1.41                 |
| R94       | 14.4       | 14.7       | 14.55             | 22.5                       | 1.94                 | 2.13                 | NA         | NA         | NA                | 10.5                       | 2.18                 | 1.37                 |
| R99       | 8.9        | 9.14       | 9.02              | 16.3                       | 1.93                 | 1.33                 | NA         | NA         | NA                | 7.4                        | 2.2                  | 1.33                 |
| R107      | 8.72       | 8.82       | 8.77              | 10                         | 1.83                 | 1.05                 | NA         | NA         | NA                | 6.3                        | 1.88                 | 0.95                 |

Only 3 out of the 14 samples got RNA Qubit readings (R1, R39, and R79). The rest of them were too low to be read on the Qubit. Confusingly, the nanodrop RNA concentration values were pretty decent. The nanodrop is considerably less accurate than the Qubit. 

R59 had a really low 260/230 value for DNA (0.13), while R17 had a really low 260/230 value for DNA (0.1). 

#### Gel 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/gel_20230911.JPG)

Even though the Qubit said that RNA values were too low, there are pretty decent bands for all samples. This likely means that there is little RNA present, but it is of high quality? Not sure, need to discuss w/ Ariana. 

Next steps for this batch of samples: 

- I have 500 uL + beads saved for future extractions for these samples. I can try to re-extract. 
	- When re-extracting, elute in less volume (ie 50 instead of 80 uL). Alternatively, I can try to pass the same li	quid through the RNA spin column to concentrate the sample.  
- Run Tapestation on the samples that failed Qubit. 
