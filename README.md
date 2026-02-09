# Vegetation Resurvey Analysis - Capria 2002-2024

*Supplementary material for MSc Thesis*

## Description

This repository contains the complete R code and documentation for the 
vegetation resurvey analysis conducted in the Capria study area 
(Tuscan-Emilian Apennines, Italy), comparing floristic composition 
between 2002 and 2024.

**Author:** Aberto Gozzi  
**Affiliation:** Università di Bologna  
**Thesis Title:** Il bosco che avanza: resurvey a Capria tra chiusura di chioma e turnover floristico
**Year:** 2026

## Methods Summary

- **Design:** Paired plot resurvey (15 plots, 2002 vs 2024)
- **Primary inference:** Presence/absence data (Jaccard distance)
- **Multivariate analysis:** NMDS, PERMANOVA (999 permutations, blocked by plot)
- **Indicator species:** IndVal with permutation tests
- **Dendrometric analysis:** Stand metrics, diameter distributions

**Software:** R 4.5.2  
**Key packages:** vegan 2.7-2, indicspecies 1.8.0, ggplot2 4.0.1

## Repository Contents

- `code/` - R analysis scripts
  - `AppendiceA1_analisi_principale.R` - Main vegetation analysis
  - `AppendiceA2_analisi_climatica.R` - Climate analysis (Walter diagrams)
  
- `docs/` - Documentation
  - `sessionInfo.txt` - R session info for reproducibility

- `outputs/` - Example outputs (selected figures)

## Requirements
```r
# Install required packages
install.packages(c("vegan", "dplyr", "tidyr", "ggplot2", 
                   "indicspecies", "permute", "patchwork", 
                   "RColorBrewer"))
```

## Usage

1. Modify `DATA_DIR` path in scripts to point to your data location
2. Ensure input data files are available
3. Run scripts:
```r
source("code/AppendiceA1_analisi_principale.R")
source("code/AppendiceA2_analisi_climatica.R")
```

## Data Availability

Input data structure and metadata are described in the thesis Appendix B.  
Raw data available on request from the author.

## Citation

If you use this code, please cite:

> GOZZI, Alberto. (2026). Vegetation resurvey analysis - Capria 2002-2024. 
> Zenodo.https://doi.org/10.5281/zenodo.18551150

Or cite the thesis:

> GOZZI, Alberto. (2026). Il bosco che avanza: resurvey a Capria tra chiusura di chioma e turnover floristico . MSc Thesis, Università di Bologna.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

**Alberto Gozzi**  
Email: alberto.gozzi@studio.unibo.it  
GitHub: @AlbertoGozzi (https://github.com/AlbertoGozzi)

## Acknowledgments

- Thesis supervisors: Francesco Maria Sabatini, PhD
- Study area: Capria, Tuscan-Emilian Apennines
- Original 2002 data: Laura Salvatori (2002)
