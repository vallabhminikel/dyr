# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`dyr` is a prion disease research project tracking disease progression across animal models (hamsters, rats, humanized PRNP knockin mice). It manages study data, computes RT-QuIC prion titers via the Spearman-Kärber method, generates survival analyses, and produces publication-ready figures.

## Running Scripts

Scripts set their working directory explicitly via `setwd('~/d/sci/src/dyr')` or `setwd('~/src/dyr')`. Run from RStudio or the command line:

```bash
Rscript src/dyr_multipanel.R    # compute titers + generate all figures and supplementary tables
```

`dyr_multipanel.R` is the only analysis script. It reads the committed TSVs in `data/analytic/` (and `data/weight_variance_analysis/`) and writes every display item to `display_items/`. The cohort data in `data/analytic/` was originally exported from the lab's Google Sheets, but the ingestion scripts are not part of this release; the TSVs are the canonical, version-controlled inputs.

## Architecture & Data Flow

```
data/analytic/*.tsv (committed study data)
            ↓
   dyr_multipanel.R
   (Spearman-Kärber titer calc, survival analysis, figures)
            ↓
   display_items/ (figure-*.png, data-s*.tsv, supplement.xlsx)
```

`dyr_multipanel.R` is the single, self-contained analysis/figure pipeline; it reads the committed TSVs and produces every publication display item.

**Step 1 – Titer calculation**: `dyr_multipanel.R` implements `spearman_karber()` for prion titer estimation from dilution-series RT-QuIC data, plus `smooth_p()` for monotonizing positivity fractions before titer calculation.

**Step 2 – Figures/analysis**: `dyr_multipanel.R` loads the TSV data, runs survival analysis (`survival` package), and produces the multi-panel figures and supplementary tables saved to `display_items/`.

## Key Data Files

| File | Description |
|------|-------------|
| `data/analytic/master.tsv` | Animal metadata (ID, sex, DOB, inoculation/death dates) |
| `data/analytic/meta.tsv` | Cohort metadata (species, genotype, inoculum, treatment, dose) |
| `data/analytic/dyr_rtquics.tsv` | RT-QuIC assay results (dilution series, positive replicates) |
| `data/analytic/weights.tsv` | Longitudinal weight measurements |
| `data/analytic/nests.tsv` | Nest-building behavioral scores |
| `data/analytic/*_mri.tsv` | MRI-cohort metadata, weights, surgeries, and blinding (master_mri, weights_mri, surgery_mri, blind_mri, meta_mri) |

## Key R Dependencies

`tidyverse`, `lubridate`, `survival`, `janitor`, `openxlsx`, `magick`, `smoother`, `beeswarm`, `broom`, `drc`

## Genotypes & Animal Models

- **WT** – wild-type (baseline susceptibility)
- **HE** – heterozygous PRNP modification
- **KO** – PRNP knockout (negative controls)
- Species: Syrian hamsters, rats, humanized PRNP knockin mice
- Primary inoculum: 263K prion strain
