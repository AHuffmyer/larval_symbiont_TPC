---
layout: post
title: MiniPrep Plus DNA/RNA extractions - test
date: '2023-07-13'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Larvae C/D 
---

# Extractions for Mcap larvae C/D Hawaii 2023 experiment 

This protocol is based on Jill’s coral adult extraction [protocol](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2022-10-25-MiniprepPlus-DNA:RNA-extractions.md), Maggie’s coral egg/larvae extraction [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/), and Kevin’s coral larvae extraction [protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/DNA-RNA-Extractions-on-P.-astreoides-larvae-BEAD-BEATING/). This post details the info about the extraction steps for the test samples that were collected in HI 2023 for A. Huffmyer's experiment on symbiont identity on *Montipora capitata* coral larval thermal performance. The github for that project is linked [here](https://github.com/AHuffmyer/larval_symbiont_TPC). 

### Samples 

These runs were only done with test samples and were done over 2 days. On 7/12/23, I extracted DNA and RNA from E7, E8, and E9. On 7/14/23, I extracted DNA and RNA from E10, E11, and E12. These samples had 100 uL of concentrated larvae and were preserved with 500 uL of DNA/RNA shield. 

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
	- Several 1000 uL, 200 uL, and 20 uL filtered pipette tips per sample.

### Protocol 

Protocol for 7/12/23 extractions: 

#### Sample prep

- Thaw samples in an ice bucket 
	- These samples look like they have ~500 uL in the tubes (as opposed to the 600 uL I expected from the 100 uL of larvae and the 500 of shield)

#### Homogenization

- Add glass beads to the thawed tubes. 
- Bead beat tubes on high for 2 minutes 
- Spin down in bench top centrifuge for 1 minute 

I added wayyy too many beads to the tubes. It was difficult to get all the liquid out after homogenization. Beads got caught in pipette tip when moving sample to new tube. I had to do one round of centrifugation and move to new tube to remove all of the glass beads from the sample. I ended up with ~300 uL of sample. 

#### Lysis 

- Add 30 uL of Pro-K digestion buffer (10:1 ratio of sample:digestion buffer) and 15 uL of Proteinase K (2:1 ratio of digestion buffer:Proteinase K) to each sample.
- Vortex 3x to mix. Incubate tubes at room temperature for 15 minutes. 
- Spin down tubes for 3 minutes at 9000 rpm. 
- Move 250 uL of sample into new tube. 
- Add 250 uL of lysis buffer to sample (1:1 ratio of buffer:sample). Vortex to mix and spin down. 
- Transfer sample into yellow DNA spin column and spin for 1 minute at 15000 rpm. 
- Move liquid collected in collection tube to labeled RNA tube. This is the RNA!

#### DNA extraction 

- Move yellow DNA column to new collection tube and discard old tube. 
- Add 400 uL of prep buffer. Spin for 1 minute at 15000 rpm. Discard liquid.
- Add 700 uL of wash buffer to column. Spin for 1 minute at 15000 rpm. Discard liquid. Do this twice.
- Spin dry for 2 minutes at 15000 rpm. Discard liquid.
- Move yellow DNA spin column to final gDNA 1.5 mL tube. 
- Add 50 uL of warm Tris to spin column. Let column sit for 5 minutes at room temperature. Spin for 1 minute at 15000 rpm. Do this twice. 
- Aliquot 12 uL into PCR tubes for QC. Store the rest of the gDNA at -20C. 

#### RNA extraction 

- Add 500 uL of 100% EtOH to RNA sample (1:1 ratio of ethanol:sample). Pipette up and down to mix. 
- Move sample into green RNA spin column. Spin for 1 minute at 15000 rpm. 
- Move gree RNA column to new collection tube and discard old tube. 
- Add 400 uL of wash buffer. Spin for 1 minute at 15000 rpm. Discard liquid.
- Make DNase mix in 1.5 mL tube. 
	- 75 uL of DNA digestion buffer x _____ samples = ____ buffer needed
	- 5 uL of DNase I x _____ samples = ____ DNase I needed
		- DNase I stored at -20. 
- Invert DNase mix gently to mix. 
- Add 80 uL of DNase mix to each spin column. Incubate at room temperature for 15 minutes 
- Add 400 uL of prep buffer. Spin for 1 minute at 15000 rpm. Discard liquid.
- Add 700 uL of wash buffer to column. Spin for 1 minute at 15000 rpm. Discard liquid. Do this twice.
- Spin dry for 2 minutes at 15000 rpm. Discard liquid.
- Move green RNA spin column to final RNA 1.5 mL tube. 
- Add 50 uL of DNA/RNA free water to spin column. Let column sit for 5 minutes at room temperature. Spin for 1 minute at 15000 rpm. Do this twice. 
- Aliquot 12 uL into PCR tubes for QC. Store the rest of the RNA at -80C. 

QC protocols are linked below. Here are the QC results (note - QC was done on the following day 7/13/23): 

| Sample ID | Qubit DNA1 | Qubit DNA2 | Qubit DNA average | Nanodrop DNA concentration | Nanodrop DNA 260/280 | Nanodrop DNA 260/230 | Qubit RNA1 | Qubit RNA2 | Qubit RNA average | Nanodrop RNA concentration | Nanodrop RNA 260/280 | Nanodrop RNA 260/230 |
| --------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- |
| E7        | 14.8       | 14.4       | 14.6              | 18.7                       | 1.95                 | 1.13                 | LOW        | LOW        | NA                | 0.6                        | \-2.13               | 0.14                 |
| E8        | 11.2       | 11.1       | 11.15             | 15.3                       | 2.07                 | 1.23                 | LOW        | LOW        | NA                | 0.4                        | \-0.69               | 1.02                 |
| E9        | 19.2       | 18.8       | 19                | 32.1                       | 1.86                 | 0.9                  | LOW        | LOW        | NA                | 1.6                        | 4.2                  | 0.2                  |

Gel image: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/gel_20230713.JPG)

Got no RNA and RNA values are too low to be detected :'(. DNA looks pretty good though. I think that the bead beating with so many beads was an issue. Tomorrow I will try with less beads. 

Protocol for 7/14/23 extractions:  

#### Sample prep

- Thaw samples in an ice bucket. Similarly to the previous extraction, these samples look like they have ~500 uL in the tubes (as opposed to the 600 uL I expected from the 100 uL of larvae and the 500 of shield).
- Add 200 uL of DNA/RNA shield to samples once they are thawed. 
- Vortex and spin down. 

#### Homogenization

- Add glass beads to the thawed tubes. Added much less than previous extraction.
- Bead beat tubes on high for 2 minutes 
- Spin down in bench top centrifuge for 1 minute 

Sample E10 does not look super homogenized, looks like there is either debris or larvae floating around still. 

#### Lysis 

- Move 400 uL to new tube, being careful to not move any glass beads with the sample. 
- Add 40 uL of Pro-K digestion buffer (10:1 ratio of sample:digestion buffer) and 20 uL of Proteinase K (2:1 ratio of digestion buffer:Proteinase K) to each sample.
- Vortex 3x to mix. Incubate tubes at room temperature for 15 minutes. 
- Spin down tubes for 3 minutes at 9000 rpm. 
- Transfer 450 uL of sample into yellow DNA spin column and spin for 1 minute at 15000 rpm. 
- Move liquid collected in collection tube to labeled RNA tube. This is the RNA!

#### DNA extraction 

Same as above. Only difference is that instead of 2 washes with 700 uL of wash buffer, I did the first wash with 700 uL and the second wash with 400 uL of wash buffer. 

#### RNA extraction 

Same as above. There was ~700 uL of RNA in the RNA tube, so I added 700 uL of 100% EtOH. Then I moved that mixture into the yellow RNA spin column and proceded with the protocol. The other difference is that instead of 2 washes with 700 uL of wash buffer, I did the first wash with 700 uL and the second wash with 400 uL of wash buffer. 

QC protocols are linked below. Here are the QC results (note - QC was done on the same day 7/14/23): 

| Sample ID | Qubit DNA1 | Qubit DNA2 | Qubit DNA average | Nanodrop DNA concentration | Nanodrop DNA 260/280 | Nanodrop DNA 260/230 | Qubit RNA1 | Qubit RNA2 | Qubit RNA average | Nanodrop RNA concentration | Nanodrop RNA 260/280 | Nanodrop RNA 260/230 |
| --------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- |
| E10       | 12.3       | 12.3       | 12.3              | 17.2                       | 2.04                 | 1.3                  | LOW        | LOW        | NA                | 3.7                        | 2.18                 | 2.7                  |
| E11       | 14         | 13.9       | 13.95             | 19.5                       | 2.04                 | 2.9                  | LOW        | LOW        | NA                | 2.7                        | 2.52                 | 3.69                 |
| E12       | 15.1       | 15         | 15.05             | 18.7                       | 2.06                 | 4.56                 | 13         | 13         | 13                | 4.6                        | 2.43                 | 2.57                 |

Gel image: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/gel_20230714.jpg)

Once again, RNA did not show up on the gel and was too low to be detected for 2 of the 3 samples. One of the samples (E12) had relatively high Qubit and Nanodrop values, but didn't show up on the gel...

Next time, I am going to follow Maggie's [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/) exactly and see if I get better RNA yields. 