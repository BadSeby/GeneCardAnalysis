---
title: 'GeneCardAnalysis: A Shiny Application for Reproducible RT-qPCR Data Analysis and Visualisation'
tags:
  - R
  - Shiny
  - RT-qPCR
  - gene expression
  - bioinformatics
  - heatmap
  - fold change
authors:
  - name: Sebastiano Di Bella
    orcid: 0000-0001-7161-303X
    affiliation: 1
  - name: Francesco Orilio
    affiliation: 1

affiliations:
  - name: Department of Precision Medicine in Medical, Surgical and Critical Care (Me.Pre.C.C.),
           University of Palermo, Palermo, Italy
    index: 1
date: 21 June 2026
bibliography: paper.bib
---

# Summary

Reverse transcription quantitative PCR (RT-qPCR) remains the gold standard
for quantifying gene expression in translational and basic research settings,
owing to its sensitivity, specificity, and low input requirements
[@livak2001; @pfaffl2001]. Despite widespread adoption, the computational
workflow from raw cycle-threshold (CT) values to publication-ready figures is
often implemented ad hoc in spreadsheets, introducing inconsistencies and
hampering reproducibility.

**GeneCardAnalysis** is an open-source R Shiny application that formalises
this workflow into a single, interactive, and reproducible tool. Starting from
a standard CT export file (comma- or semicolon-separated), the app performs
delta-CT normalisation [@livak2001], computes log2 fold-change (log2FC)
relative to a user-specified reference condition, and calculates per-gene
z-scores on the log2 scale. Results are rendered as publication-quality
heatmaps via ComplexHeatmap [@gu2016] and, when only two conditions are
compared, as ranked bar charts. All outputs — PNG, PDF, Excel, and ZIP
archives — can be downloaded directly from the interface.

Experimental comparisons are defined in a YAML configuration file that
travels alongside the data, making the analysis fully reproducible and
shareable. An in-app form allows users without YAML experience to configure
comparisons interactively and export the resulting YAML for future reuse.


# Statement of Need

Wet-lab researchers routinely generate RT-qPCR datasets containing dozens of
genes across multiple experimental conditions. The standard delta-CT analysis
[@livak2001] is conceptually straightforward, but its implementation is
error-prone when done manually: incorrect reference gene averaging, accidental
z-scoring of linear (2^(-deltaCT)) rather than log2 values, and inconsistent
handling of "Undetermined" calls are common pitfalls.

Dedicated commercial software (e.g., QuantStudio Design & Analysis Software,
Bio-Rad CFX Manager) automates part of this workflow but produces proprietary
outputs that are difficult to integrate with downstream R-based analyses.
Existing open-source tools such as `HTqPCR` [@perkins2012] and `pcr`
[@ahmed2018] are command-line packages that require R programming fluency and
do not provide interactive visualisation.

GeneCardAnalysis fills this gap by offering:

1. A **no-code interface** accessible to biologists without programming experience.
2. **Methodologically correct** z-scoring on the log2 scale (= -deltaCT), not
   on the exponential 2^(-deltaCT) values.
3. **YAML-driven reproducibility** — the full analysis specification is
   captured in a human-readable file that can be version-controlled.
4. **Publication-ready outputs** generated with a single click, including
   high-resolution PNG and vector PDF heatmaps annotated with project and
   author metadata.


# Functionality

## Input

The application accepts a CSV file with at minimum three columns:
`Sample.Name`, `Target.Name`, and `CT`. This format is compatible with direct
exports from Applied Biosystems QuantStudio and Bio-Rad CFX systems. The
separator (comma or semicolon) is auto-detected. "Undetermined" calls are
replaced with a user-configurable numeric value (default: 40).

## Analysis pipeline

1. **Housekeeping normalisation.** CT values are averaged across technical
   replicates and corrected by subtracting the mean CT of the designated
   housekeeping gene(s) (GAPDH and/or HPRT1, or user-specified). The result
   is converted to a normalised linear expression value: 2^(-deltaCT).

2. **Log2 fold-change.** For each comparison, log2FC is computed as
   log2(sample / reference), where the reference is the user-specified
   control condition (e.g., untreated, UT).

3. **Z-score.** Per-gene z-scores are calculated across samples on the log2
   scale (log2(2^(-deltaCT)) = -deltaCT), preserving the symmetric, additive
   properties of log-transformed data.

## Visualisation

- **Z-score heatmap** — diverging colour scale centred at zero; suitable for
  comparing relative expression patterns across many samples.
- **Normalised heatmap** — colour scale based on empirical quantiles of
  2^(-deltaCT) values; suitable for comparing absolute expression levels.
- **Log2FC heatmap / bar chart** — when three or more samples are present, a
  heatmap diverging around zero is shown; with two conditions, the tool
  switches automatically to a ranked horizontal bar chart for clarity.

All three visualisations support user-defined colour palettes, configurable
row/column clustering (none, rows, columns, or both), and manual column
ordering.

## Configuration and reproducibility

Comparisons are specified in a YAML file with the following structure:

```yaml
project_name: EMT_study
author: SDB
global_settings:
  housekeeping_genes: [GAPDH, HPRT1]
  undetermined_value: 40
  exclude_genes: [RT, gDNA]
comparisons:
  - name: EMT_vs_UT
    reference: UT
    samples: [EMT_24h, EMT_48h, UT]
    signature_type: signature_emt
    cluster_heatmap: true
```

Custom gene signatures can be embedded in the same YAML under
`custom_signatures`. The full configuration can also be generated via the
in-app form and exported as YAML.

## Gene signatures

Four built-in gene panels are included (all genes, EMT, metastasis, and
metabolism), derived from curated literature panels used in the Stassi Lab.
Custom signatures defined in the YAML are loaded dynamically and appear
alongside built-in options in the interface.


# Availability

GeneCardAnalysis is implemented in R [@rcoreteam] and is available at
<https://github.com/BadSeby/GeneCardAnalysis> under the MIT License.
The application can be run locally with:

```r
# Install dependencies
install.packages(c("shiny", "bslib", "circlize", "readxl", "tidyverse",
                   "grid", "glue", "writexl", "showtext", "yaml",
                   "DT", "colourpicker", "gtools"))
BiocManager::install("ComplexHeatmap")

# Launch
shiny::runApp("app.R")
```

A hosted demo instance is available at: <!-- TODO: add shinyapps.io URL -->


# Acknowledgements

The authors thank the members of the Stassi Lab for testing the application
and providing feedback on clinical workflows. <!-- TODO: add funding statement -->


# References
