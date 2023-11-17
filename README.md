[![DOI: 10.1093/bib/bbaa269](https://img.shields.io/badge/DOI-10.1101/2023.09.08.556895-brightgreen)](https://doi.org/10.1101/2023.09.08.556895)
# MUSTANG
![GitHub Logo](/Miscel/Fig1_A4_cropped.png)
## Overview
**MU**lti-sample **S**patial **T**ranscriptomics data **AN**alysis with cross-sample transcriptional similarity **G**uidance (**MUSTANG**) is a computaional framework, which is capable of performing multi-sample spatial transcriptomics spot cellualar deconvolution by allowing both cross-sample expression based similarity information sharing as well as spatial correlation in gene expression patterns within samples.

## Citation

If you use this code, please cite the [preprint](https://www.biorxiv.org/content/10.1101/2023.09.08.556895v1):

```
MUSTANG: MUlti-sample Spatial Transcriptomics data ANalysis with cross-sample transcriptional similarity Guidance
Seyednami Niyakan, Jianting Sheng, Yuliang Cao, Xiang Zhang, Zhan Xu, Ling Wu, Stephen T.C. Wong, Xiaoning Qian
bioRxiv 2023.09.08.556895; doi: https://doi.org/10.1101/2023.09.08.556895
```
## Quick Start
In order to analyze your multi-sample spatial transcriptional (ST) data with MUSTANG, 4 main steps should be performed:

1.  **Spots Spatial Graph**: The adjacency matrix of spots spatial graph should be extracted based on the layout.
1.  **Spots Transcriptional Graph**: The adjacency matrix of spots transcriptional graph in which spots that are transcriptionally similar to eachother are connected with an edge should be extracted.
1.  **Spots Similarity Graph**: The adjacency matrix of spots similarity graph should be constructed based on adjacency matrices of spots spatial and transcriptional graphs. 
1.  **Bayesian Deconvolution Analysis**: The Poisson discrete deconvolution model should be applied to extract the deconvolution parameters.

## Tutorials
- [Analysis of Mouse Brain ST data with `MUSTANG`](https://github.com/namini94/MUSTANG/blob/main/Tutorial/Mouse%20Brain%20/Mouse_Brain.md)

