

**PrP lowering is effective against prion disease across myriad in vivo models**

This repo contains analysis code, source data, and display items for a manuscript in preparation.

## Overview

Here you will find everything needed to reproduce the figures and supplementary
data tables from the manuscript. It tracks prion disease progression across multiple
animal models — prion-inoculated rats and Syrian hamsters, spontaneously sick mice
overexpressing bank vole PrP (Tg(Vole)), and humanized mice infected
with human prion isolates — and quantifies the effect of lowering PrP, either
genetically or with antisense oligonucleotides (ASOs), on survival, body weight,
plasma neurofilament light (NfL), brain/CSF prion seeding activity (RT-QuIC).

## Repository layout

```
src/
  dyr_multipanel.R                 # the single analysis + figure pipeline
data/
  analytic/                        # cleaned, version-controlled inputs (TSV)
  weight_variance_analysis/        # baseline weight-volatility sub-analysis + raw data
  images/                          # histology / reference images used in figures
  sequencing/                      # knockout-allele genotyping reads and PrP reference sequences
display_items/
  figure-1.png … figure-4.png      # main figures
  figure-s1.png … figure-s5.png    # supplementary figures
  data-s001.tsv, data-s002.tsv, …  # supplementary / source data (individual files)
  supplement.xlsx                  # supplementary / source data (combined workbook)
CLAUDE.md                          # developer notes on architecture and data flow
```

## Reproducing the analysis

The entire analysis is a single, self-contained R script. From the repository root:

```bash
Rscript src/dyr_multipanel.R
```

It reads the committed TSVs in `data/analytic/` (and
`data/weight_variance_analysis/`), computes RT-QuIC prion titers by the
Spearman-Kärber method, runs the survival and biomarker analyses, and writes out
figures (`display_items/figure-*.png`), supplementary datasets
(`display_items/data-s*.tsv`), and the combined workbook
(`display_items/supplement.xlsx`).

The script sets its working directory explicitly (`setwd('~/d/sci/src/dyr')` /
`setwd('~/src/dyr')`); if running in interactive mode, then edit this to match your local path before running. In batch mode this is not relevant.

### R dependencies

`tidyverse`, `lubridate`, `survival`, `janitor`, `openxlsx`, `magick`, `smoother`,
`beeswarm`, `broom`, `drc`

```r
install.packages(c("tidyverse","lubridate","survival","janitor","openxlsx",
                   "magick","smoother","beeswarm","broom","drc"))
```

## Data

The cohort data in `data/analytic/` was originally exported from the study's
electronic records; those exported TSVs are the canonical, version-controlled inputs
for this repository, and no live database or network access is required to reproduce
the results.

Genotype, dose, and treatment definitions for each cohort are in
`data/analytic/meta.tsv` and `data/analytic/cohort_lookup.tsv`. Per-animal records are
in `data/analytic/master.tsv` (and `master_mri.tsv` for the MRI-facility cohorts).

## Citation

This is a manuscript in progress, no citation is available yet

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.


## Contact

<!-- TODO: Add corresponding-author contact / lab info. -->
