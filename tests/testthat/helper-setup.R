# ==============================================================================
# Shared test helpers — loaded automatically by testthat before any test file
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(gtools)
  library(ggplot2)
})

# During R CMD check the package is installed and functions are in the namespace.
# For local interactive testing (testthat::test_dir), source from R/ if needed.
if (!exists("extract.normalised.CTs", mode = "function")) {
  pkg_root <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = FALSE)
  r_dir    <- file.path(pkg_root, "R")
  if (dir.exists(r_dir)) {
    source(file.path(r_dir, "card_functions.R"),        local = FALSE)
    source(file.path(r_dir, "signatures_genes_card.R"), local = FALSE)
  }
}

# ------------------------------------------------------------------------------
# Toy CT data builders
# ------------------------------------------------------------------------------

#' Minimal 2-sample, 3-gene dataset (single replicates)
#' UT:  GAPDH=20, GENE1=25 (deltaCT=5 → 2^-5=0.03125),
#'               GENE2=30 (deltaCT=10 → 2^-10≈9.766e-4)
#' EMT: GAPDH=20, GENE1=23 (deltaCT=3 → 2^-3=0.125),
#'               GENE2=28 (deltaCT=8  → 2^-8≈3.906e-3)
make_cts_simple <- function() {
  data.frame(
    Sample.Name = rep(c("UT", "EMT"), each = 3),
    Target.Name = rep(c("GAPDH", "GENE1", "GENE2"), times = 2),
    CT          = c(20, 25, 30,   20, 23, 28),
    stringsAsFactors = FALSE
  )
}

#' Dataset with comma decimal separators in CT column
make_cts_comma <- function() {
  df <- make_cts_simple()
  df$CT <- gsub("\\.", ",", as.character(df$CT))
  df
}

#' Dataset with an "Undetermined" call in GENE2/UT
make_cts_undetermined <- function() {
  df <- make_cts_simple()
  df$CT[df$Sample.Name == "UT" & df$Target.Name == "GENE2"] <- "Undetermined"
  df$CT <- as.character(df$CT)
  df
}

#' Dataset with technical duplicates (all genes appear twice per sample)
make_cts_duplicates <- function() {
  df <- make_cts_simple()
  # Add jitter to second replicate
  df2 <- df
  df2$CT <- df2$CT + c(0.1, 0.2, 0.3,  0.1, 0.2, 0.3)
  rbind(df, df2)
}

#' Dataset with partial duplicates (only GENE1 is duplicated per sample)
make_cts_partial_dup <- function() {
  df <- make_cts_simple()
  extra <- df[df$Target.Name == "GENE1", ]
  extra$CT <- extra$CT + 0.5
  rbind(df, extra)
}

#' Dataset with both GAPDH and HPRT1 as housekeeping genes
make_cts_two_hk <- function() {
  base <- data.frame(
    Sample.Name = rep(c("UT", "EMT"), each = 4),
    Target.Name = rep(c("GAPDH", "HPRT1", "GENE1", "GENE2"), times = 2),
    CT          = c(20, 21, 25, 30,   20, 21, 23, 28),
    stringsAsFactors = FALSE
  )
  base
}

#' Dataset with an exclude-control gene (gDNA) present
make_cts_with_control <- function() {
  df <- make_cts_simple()
  ctrl <- data.frame(
    Sample.Name = c("UT", "EMT"),
    Target.Name = c("gDNA", "gDNA"),
    CT          = c(38, 39),
    stringsAsFactors = FALSE
  )
  rbind(df, ctrl)
}

#' Dataset where housekeeping gene is completely absent
make_cts_no_hk <- function() {
  data.frame(
    Sample.Name = c("UT", "UT"),
    Target.Name = c("GENE1", "GENE2"),
    CT          = c(25, 30),
    stringsAsFactors = FALSE
  )
}

#' Dataset that triggers auto-select: HPRT1 has lower variance than GAPDH
make_cts_autoselect_hprt1 <- function() {
  data.frame(
    Sample.Name = rep(c("S1", "S2", "S3"), each = 3),
    Target.Name = rep(c("GAPDH", "HPRT1", "GENE1"), times = 3),
    CT          = c(
      20, 21, 28,   # S1: GAPDH var=high across samples; HPRT1 var=low
      22, 21, 26,   # S2
      24, 21, 27    # S3
    ),
    stringsAsFactors = FALSE
  )
}
