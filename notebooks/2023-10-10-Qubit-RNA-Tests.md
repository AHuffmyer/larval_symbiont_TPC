---
layout: post
title: RNA Qubit tests
date: '2023-10-10'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Larvae C/D 
---

### Qubit RNA tests 

Since I've been getting NAs for a lot of the RNA Qubit values, I'm going to try a few different iterations of Qubit today: 

- 1) BR Qubit w/ 1 uL sample input
- 2) BR Qubit w/ 2 uL sample input
- 3) BR Qubit w/ 1 uL sample input - Prada kit 
- 4) HS Qubit w/ 1 uL sample input
- 5) HS Qubit w/ 2 uL sample input

I'm going to use 6 samples: 

- R5 (a sample that passed Qubit but has not been run on the Tapestation; eluted in 80 uL)
- R20 (a sample that passed Qubit but has not been run on the Tapestation; eluted in 65 uL)
- R43 (a sample that has failed Qubit but had high concentration on the Tapestation; eluted in 65 uL)
- R54 (a sample that has failed Qubit but had high concentration on the Tapestation; eluted in 80 uL)
- R77 (a sample that passed Qubit but has not been run on the Tapestation; eluted in 65 uL)
- R107 (a sample that has failed Qubit but had low concentration on the Tapestation; eluted in 80 uL)


R5 was extracted on [7/25/23](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-07-25-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae.md). R54 was extracted on [9/16/23](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-09-16-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae.md) and Tapestation was run on [10/5/23](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-10-05-Tapestation-McapLarvae.md). R107 was extracted on [9/11/23](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-09-11-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae.md) and Tapestation was run on [9/22/23](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-09-22-Tapestation-McapLarvae.md). These samples were eluted in 80 uL. 

R20, R43, and R77 were extracted on [10/5/23](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-10-05-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae.md). R43 was run on the Tapestation on [10/5/23](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-10-05-Tapestation-McapLarvae.md). These samples were eluted in 65 uL. 

**Important note:** If the sample input is 2 uL, the input volume must be increased on the Qubit from 1 ng/uL to 2 ng/uL before running. 

#### Results 

The plot unfortunately thickens. 

| Sample ID | Elution vol (uL) | Qubit BR [ ] | Nanodrop [ ] | Tapestation [ ] | Qubit BR [ ] 1uL | Qubit BR [ ] 2uL | Qubit BR [ ] 1uL PRADA | Qubit HS [ ] 1uL | Qubit HS [ ] 2uL |
| --------- | ---------------- | ------------ | ------------ | --------------- | ---------------- | ---------------- | ---------------------- | ---------------- | ---------------- |
| R5        | 80               | 19.2         | 9            | Not run         | NA               | 6.95             | 14.2                   | 8.11             | 7.14             |
| R54       | 80               | NA           | 7.8          | 18              | NA               | 7.4              | 12.1                   | 7.95             | 7.4              |
| R107      | 80               | NA           | 6.3          | 8.13            | NA               | 5.45             | NA                     | 5.43             | 5.6              |
| R20       | 65               | 13           | 11.3         | Not run         | 12               | 11.2             | 17.3                   | 10.02            | 10.05            |
| R43       | 65               | NA           | 8.5          | 13.3            | NA               | 8.15             | 13.8                   | 6.99             | 7.5              |
| R77       | 65               | 11.8         | 8.9          | Not run         | NA               | 8.55             | 16                     | 8.93             | 7.7              |

Columns: 

- Sample ID 
- Elution vol (uL) 
- Qubit BR - original Qubit concentration from the day that sample was extracted
- Nanodrop - original Nanodrop concentration from the day that sample was extracted
- Tapestation - concentration from TS run 
- BR Qubit 1 uL - BR Qubit w/ 1 uL sample input
- BR Qubit 2 uL - BR Qubit w/ 2 uL sample input
- BR Qubit 1 uL PRADA - BR Qubit w/ 1 uL sample input using the Prada Qubit kit 
- HS Qubit 1 uL - HS Qubit w/ 1 uL sample input
- HS Qubit 2 uL - HS Qubit w/ 2 uL sample input

#### Thoughts

Even though R5 and R77 got relatively high concentrations from the original Qubit value, both samples got an NA when re-ran with the BR kit with a sample input of 1 uL. I got much higher BR values using the Prada lab kit, with the exception of R107. I got similar values for the HS kit with sample inputs of 1 uL & 2 uL. The Tapestation values for R54 is very high compared to the Qubit concentrations from both the BR and HS kits. 

The Tapestation concentration for R43 is similar to the Qubit concentration obtained by the Prada lab kit. 

#### Next steps 

- Continue eluting in 65 uL
- Run Qubit BR using 2 uL of sample input
- Discuss w/ Ariana and Hollie 