---
layout: post
title: ITS2 clean-up test for Mcap C/D
date: '2024-04-17'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Larvae C/D 
---

# ITS2 for Mcap larvae C/D Hawaii 2023 experiment 

This post details information on ITS2 sequencing for the Mcap larval C/D 2023 experiment. The github for that project is linked [here](https://github.com/AHuffmyer/larval_symbiont_TPC).  ITS2 amplicon sequencing is being used to characterize the symbiont community from larval samples. See Ariana's [post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/ITS2-amplicon-PCR-and-preparation-for-sequencing-20240326/) about the ITS2 protocol for this project. 

### Samples 

Samples (18 samples, positive control, negative control) from A. Huffmyer's 2023 experiment were amplified but there appeared to be some contamination in the gels, as evidenced by bands ~100bp. After further troubleshooting, it was determined that these bands are likely [primer dimers](https://www.minipcr.com/primer-dimer-pcr/). Therefore, we decided to do a bead clean-up to remove the primer dimer in the samples. Today, I ran a test clean-up. Samples used today were: 

- R58
- R65
- R72
- Positive control (Mcap 2020 sample)
- Negative control (ultrapure h2o

Samples were thawed on ice. When they thawed, there was a weird-looking separation happening in all of the tubes (see photo below). I'm not sure what the bottom layer was, but I assumed that the PCR product was the top layer. Flicking and brief vortexing did not seem to affect it. For the protocol, I took material from the top layer.

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/samples_20240417.jpg)

### Equipment and materials 

- [Kapa Pure Beads](https://elabdoc-prod.roche.com/eLD/web/pi/en/products/SEQ-KAPA-0161?searchTerm=07983271001&catalog=Researcher&orderBy=Relevance)
- Tris-HCl (10 mM)
- Molecular grade ethanol
- Magnetic stand 
- PCR tubes 
- Gel rig + materials 

### Protocol 

Protocol was followed according to this [post](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/_posts/2024-04-01-KAPA-bead-clean-up-protocol-for-removal-of-primer-dimers-from-PCR-product.md). A few notes: 

- After adding Kapa Beads, I incubated tubes for 15 minutes at room temperature.
- After the ethanol washes, I dried beads for ~4 minutes or until beads were matte/ethanol had evaporated. 
- After adding the elution buffer (10 mM Tris-HCl), I incubated tubes for 8 minutes at room temperature. 

### QC 

Following the clean-up protocol, I ran a 2% gel for 90 minutes at 80 volts to check if the primer dimers were removed. I also used a 1kb ladder and a 100 bp ladder. 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_20240417.JPG)

There still looks to be some faint primer dimer contamination. Maybe another bead clean-up or more ethanol wash steps would help? Will send to Hollie and Ariana for their thoughts. 

