# shotgun-patchSCP
# Patch-SCP: Linking Electrophysiology and Proteomics in Brain Slices

This repository contains the analysis code used in  
**"Patch-Clamp Single-Cell Proteomics in Acute Brain Slices: A Framework for Recording, Retrieval, and Interpretation" (Biorxiv DOI: https://doi.org/10.1101/2025.09.15.675920, 2025, under review).**

---

## Contents
- Supplementary CSV/XLSX tables (Zip file)
- SynGO Analyses (Zip file)
- Example output plots (Zip file)
- R scripts for figure generation (Figures 2–6 and S2-3)
  
| Figure (paper) | Script file name |
|----------------|------------------|
| Figure 2 – Protein overlap (Venn) | `Figure 2 Venn Diagram.R` |
| Figure 3B – Gigaseal preserved (BP Upset) | `Figure 3B Gigaseal preserved upset for BPs.R` |
| Figure 3C – Gigaseal preserved (BP Heatmap) | `Figure 3C Gigaseal preserved BP Heatmap.R` |
| Figure 5 – Split SynGO heatmaps | `Figure 5 Split SynGO heatmap 9-23-25.R` |
| Figure 6 – Ion channel heatmap (toggle for duplicates) | `Figure 6 Channel Heatmap with toggle 9-26-25.R` |
| Supplementary Figure S2A – BP Upset | `S2A BP Upset 9-26-25.R` |
| Supplementary Figure S2B – CC Upset | `S2B CC Upset 9-26-25.R` |
| Supplementary Figure S3 – GPCR heatmap (toggle) | `S3 GPCR Heatmap with toggle 9-26-25.R` |

## Data Availability
Raw and processed mass spectrometry data can be found using the MassIVE dataset identifier MSV000099156 and ProteomeXchange dataset identifier PXD068359. Supplementary tables S1–S3 are provided with the manuscript.

## Usage
Note that all scripts require report.pg_matrix.tsv and SynGO analysis files. If raw spectra files and/or SynGO outputs are being analyzed differently, modify names in the metadata file. Then open R and run the scripts directly (e.g., `source("Figure2_Venn_Diagram.R")`). Scripts are annotated with input file paths and parameters. Script for Figure 5 must be run first. 

## Citation
If you use these scripts, please cite:
Rodriguez et al., *Biorxiv*, 2025 (DOI: https://doi.org/10.1101/2025.09.15.675920).

## Dependencies
All analyses were performed in **R Version 4.3.3** with the following packages:
- `readr`
- `readxl`
- `dplyr`
- `stringr`
- `tidyr`
- `lubridate`
- `ggplot2`
- `eulerr`
- `pheatmap`
- `ComplexHeatmap`
- `circlize`
- `RColorBrewer`
  
To install all dependencies in one step:  

```r
install.packages(c("readr","readxl","dplyr","stringr","tidyr",
                   "lubridate","ggplot2","eulerr","pheatmap","RColorBrewer"))
BiocManager::install(c("ComplexHeatmap","circlize"))

