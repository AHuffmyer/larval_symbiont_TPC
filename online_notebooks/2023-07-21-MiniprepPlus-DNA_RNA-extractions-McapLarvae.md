---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2023-07-21'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Larvae C/D 
---

# Extractions for Mcap larvae C/D Hawaii 2023 experiment 

This protocol is based on Maggie’s coral egg/larvae extraction [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/) and Kevin’s coral larvae extraction [protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/DNA-RNA-Extractions-on-P.-astreoides-larvae-BEAD-BEATING/). This post details the info about the extraction steps for the test samples that were collected in HI 2023 for A. Huffmyer's experiment on symbiont identity on *Montipora capitata* coral larval thermal performance. The github for that project is linked [here](https://github.com/AHuffmyer/larval_symbiont_TPC). 

### Samples 

The samples run today were R54, R81, and R104. These samples had 50 larvae in them and were preserved with 700 uL of DNA/RNA Shield. This photo is post-bead beating: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/samples_20230721.JPG)

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

#### Sample prep

- Thaw samples in an ice bucket. These samples have 50 larvae in 700 uL of DNA/RNA shield. 

#### Homogenization

- Spin down tubes so there are no bubbles in lid 

- Add 0.25 uL of glass beads (0.5mm) into sample tubes. 

- Beat beat tubes at high speed for 2 minutes 

- Spin down in mini centrifuge for 1 minute 

#### Lysis 

- Remove 500 uL from sample and add it to new tube. Save the leftover sample + beads at -80°C. 

- Add 50 uL of Proteinase K digestion buffer (10:1 ratio of sample:digestion buffer) and 25 uL of Proteinase K (2:1 ratio of digestion buffer:Proteinase K) to each sample

- Vortex and spin down. 

- Add 575 uL of lysis buffer (1:1 ratio of sample:lysis buffer) to each sample. The volume in each tube should be 1150 uL. 

- Vortex and spin down.  

- Transfer 700 uL of the sample into the corresponding yellow DNA spin column. Spin for 30 seconds at 16000 rpm.

- Move the liquid collected in the collection tube into a labeled 5 mL tube. This is your RNA!  

- Transfer the remainder of the sample (~450 uL) into the yellow DNA spin column. Spin for 30 seconds at 16000 rpm.

- Move the liquid collected in the collection tube into a labeled 5 mL tube. This is your RNA! Set this aside until RNA extraction steps. 

#### DNA extraction

- Move yellow DNA column to new collection tube and discard the old one. 

- Add 400 uL of prep buffer. Spin for 30 seconds at 16000 rpm. Discard liquid. 

- Add 700 uL of wash buffer. Spin for 30 seconds at 16000 rpm. Discard liquid. 

- Add 400 uL of wash buffer. Spin for 2 minutes at 16000 rpm. Discard liquid. 

- Spin columns dry for 2 minutes at 16000 rpm. 

- Move yellow DNA spin column into final gDNA 1.5 mL tube. 

- Add 50 µL of warm Tris to spin column. Let columns sit for 5 minutes at room temperature. Spin for 30 seconds at 16000 rpm. DO THIS TWICE!

- Aliquot 12 µL into PCR tubes for QC. Store the rest of the gDNA at -20℃. 

#### RNA extraction

- Add 1150 uL of 100% EtOH (1:1 ratio of sample:ethanol) to the 5 mL tube with the RNA. Pipette up and down to mix. 

- Transfer 700 uL of sample into the green RNA spin column. Spin for 30 seconds at 16000 rpm. Discard liquid. 

- Repeat the above step until all liquid has been moved from the 5 mL tube. 

- Add 400 uL of wash buffer. Spin for 30 seconds at 16000 rpm. Discard liquid. 

- Make DNase rxn mix 
	- To make DNase rxn mix: 
		- DNA digestion buffer: 75 µL x _____ samples = _____ buffer needed
		- DNase I: 		 5 µL x _____ samples = _____ DNase I needed
	- DNase I is stored at -20℃
	- Combine in 1.5 mL tube and invert gently to mix

- Add 80 uL of DNase mix to each spin column. Incubate at room temperature for 15 minutes. 

- Add 400 uL of prep buffer. Spin for 30 seconds at 16000 rpm. Discard liquid. 

- Add 700 uL of wash buffer. Spin for 30 seconds at 16000 rpm. Discard liquid. 

- Add 400 uL of wash buffer. Spin for 2 minutes at 16000 rpm. Discard liquid. 

- Spin columns dry for 2 minutes at 16000 rpm. 

- Move green RNA spin column into final RNA 1.5 mL tube. 

- Add 50 µL of warmed DNA/RNA free water to spin column. Let columns sit for 5 minutes. Spin for 30 seconds at 16000 rpm. DO THIS TWICE!

- Aliquot 12 µL into PCR tubes for QC. Store the rest of the RNA at -80℃. 

### QC 

#### Qubit & Nanodrop results 

| Sample ID | Qubit DNA1 | Qubit DNA2 | Qubit DNA average | Nanodrop DNA concentration | Nanodrop DNA 260/280 | Nanodrop DNA 260/230 | Qubit RNA1 | Qubit RNA2 | Qubit RNA average | Nanodrop RNA concentration | Nanodrop RNA 260/280 | Nanodrop RNA 260/230 |
| --------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- |
| R56       | 24.2       | 24         | 24.1              | 47.2                       | 1.99                 | 1.98                 | 19.4       | 19.2       | 19.3              | 12.1                       | 1.88                 | 0.95                 |
| R81       | 18.1       | 18.1       | 18.1              | 34.6                       | 1.69                 | 1.69                 | 15.2       | 14.8       | 15                | 9.5                        | 1.95                 | 0.81                 |
| R104      | 20.6       | 20.6       | 20.6              | 34.2                       | 2.02                 | 2.11                 | 15.6       | 15.4       | 15.5              | 8.5                        | 2.03                 | 0.66                 |

Finally got some good RNA concentrations!! Woohoo!

#### Gel part 1

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/gel_20230721.JPG)

Gel has 2 bands in RNA side (right side). I can also see 2 bands that look like RNA on the DNA side (left side). Did I accidently contaminate the DNA w/ RNA? Hard to imagine that because I extracted the DNA first and then aliquoted it out before starting RNA extractions. 

#### Gel part 2 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/gel_20230724.JPG)

This gel was done on 7/24/23, but it looks similar to the gel above.

#### Tapestation 

Tapestation run on 7/24/23. 

##### DNA 

- Link to PDF of DNA results [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_Mcap2023_2023-07-24.pdf)

Overview of DNA results: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20230724.png)

R54 DNA results: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_R54_20230724.png)

R81 DNA results: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_R81_20230724.png)

R104 DNA results: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_R104_20230724.png)

Overall, DNA looks good. One band at the top of the TS gel, which is to be expected. 

##### RNA 

- Link to PDF of RNA results [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_Mcap2023_2023-07-24.pdf)

Overview of RNA results: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_overview_20230724.png)

R54 RNA results: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_R54_20230724.png)

R81 RNA results: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_R81_20230724.png)

R104 RNA results: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_R104_20230724.png)

R81 and R104 did not have a RIN score calculated because the concentration was >10 ng/ul. That's okay - just because a RIN couldn't be calculated doesn't mean there isn't RNA in there. The TS gel shows bands at 18S and 28S for all samples, which is good. There are also some bands between 25-200 nt in the TS gel, which may be due to lots of small RNAs. There are peaks in R81 and R104 that are not 18S or 28S, which may also be due to small RNAs and/or contimination. 

### Notes moving forward 

- Talked to Jacob Green (Puritz lab) while doing TS. He recommended mixing the RNA (either by vortexting, inverting, etc) to get higher RNA concentrations. To get a higher concentration of RNA, he also recommended either eluting in a lower volume (i.e., 50 uL instead of 100 uL DNA/RNA free water) or passing the same volume through the RNA column twice. 
- We still have decent RNA concentrations. I may add another wash w/ wash buffer before eluting the RNA. 

