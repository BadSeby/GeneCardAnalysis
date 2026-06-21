library(gtools)
library(ggplot2)
library(dplyr)

# ==============================================================================
# extract.normalised.CTs
# ------------------------------------------------------------------------------
# Normalises CT values from RT-qPCR data using the delta-CT method.
#
# Parameters:
#   CTs               data.frame with columns: Target.Name, Sample.Name, CT
#   undetermined      numeric value to substitute for "Undetermined" (default 40)
#   housekeeping.gene vector of housekeeping gene names (default: auto-selects
#                     GAPDH or HPRT1 based on minimum variance)
#   exclude_genes     vector of genes to exclude (controls, non-target wells).
#                     Default: c("RT", "gDNA", "GDNA", "CTRN", "RQ1", "RQ2", "PCR")
#   handle.NAs        if TRUE, replaces NAs with values near the minimum (default TRUE)
#
# Returns:
#   numeric matrix genes x samples with 2^(-deltaCT) values
# ==============================================================================
extract.normalised.CTs <- function(CTs,
                                   undetermined      = 40,
                                   housekeeping.gene = NULL,
                                   exclude_genes     = c("RT", "gDNA", "GDNA",
                                                         "CTRN", "RQ1", "RQ2", "PCR"),
                                   handle.NAs        = TRUE) {

  CTs <- CTs[, c("Target.Name", "Sample.Name", "CT")]
  CTs[which(CTs$CT == "Undetermined"), "CT"] <- undetermined
  CTs$CT <- gsub(",", ".", CTs$CT)   # handle comma as decimal separator
  CTs$CT <- as.numeric(CTs$CT)

  # Automatic housekeeping gene selection (GAPDH vs HPRT1) if not specified
  if (is.null(housekeeping.gene)) {
    gapdh_var  <- CTs %>% filter(Target.Name == "GAPDH") %>% pull(CT) %>% var(na.rm = TRUE)
    hprt1_var  <- CTs %>% filter(Target.Name == "HPRT1") %>% pull(CT) %>% var(na.rm = TRUE)
    housekeeping.gene <- ifelse(!is.na(gapdh_var) && !is.na(hprt1_var) && gapdh_var < hprt1_var,
                                "GAPDH", "HPRT1")
    message(sprintf("[Housekeeping] '%s' selected automatically (minimum variance).",
                    housekeeping.gene))
  }

  # Genes to remove: exclude_genes supplied by user,
  # MINUS the housekeeping genes (which must remain for delta-CT calculation)
  genes_to_remove <- setdiff(exclude_genes, housekeeping.gene)
  CTs <- CTs[which(!(CTs$Target.Name %in% genes_to_remove)), ]

  samples <- unique(CTs$Sample.Name)

  CT.matrix <- lapply(samples, function(sample) {
    columns <- CTs[which(CTs$Sample.Name == sample), ]

    # Average duplicates (if all genes are duplicated)
    tutti_duplicati  <- all(columns$Target.Name %in%
                            columns$Target.Name[duplicated(columns$Target.Name)])
    alcuni_duplicati <- any(columns$Target.Name %in%
                            columns$Target.Name[duplicated(columns$Target.Name)])

    if (tutti_duplicati) {
      columns <- aggregate(columns$CT,
                           by  = list(Target.Name = columns$Target.Name,
                                      Sample.Name = columns$Sample.Name),
                           FUN = mean)
      colnames(columns)[3] <- "CT"
    } else if (alcuni_duplicati) {
      warning(sprintf("[Sample: %s] Not all genes are duplicated.", sample))
      # Average whatever duplicates exist so row names remain unique
      columns <- aggregate(columns$CT,
                           by  = list(Target.Name = columns$Target.Name,
                                      Sample.Name = columns$Sample.Name),
                           FUN = mean)
      colnames(columns)[3] <- "CT"
    }

    rownames(columns) <- columns$Target.Name
    columns$Target.Name <- NULL
    colnames(columns)[which(colnames(columns) == "CT")] <- sample
    columns$Sample.Name <- NULL

    # Check housekeeping gene presence in the sample
    hk_present <- intersect(housekeeping.gene, rownames(columns))
    if (length(hk_present) == 0) {
      stop(sprintf(
        "No housekeeping gene found in sample '%s'. Expected: %s",
        sample, paste(housekeeping.gene, collapse = ", ")
      ))
    }
    if (length(hk_present) < length(housekeeping.gene)) {
      missing_hk <- setdiff(housekeeping.gene, hk_present)
      warning(sprintf("[Sample: %s] Missing housekeeping genes: %s. Used: %s.",
                      sample,
                      paste(missing_hk, collapse = ", "),
                      paste(hk_present,  collapse = ", ")))
    }

    # Delta-CT: subtract the mean of housekeeping genes
    columns[, 1] <- columns[, 1] - mean(columns[hk_present, 1], na.rm = TRUE)

    # Remove housekeeping genes from the output matrix
    columns <- columns[which(!(rownames(columns) %in% housekeeping.gene)), , drop = FALSE]

    # 2^(-deltaCT) — required for downstream heatmap analysis
    columns[, 1] <- 2^(-columns[, 1])
    return(columns)
  })

  CT.matrix <- do.call("cbind", CT.matrix)
  CT.matrix <- as.matrix(CT.matrix[order(row.names(CT.matrix)), ])

  # NA handling: replace with random values near the minimum
  if (isTRUE(handle.NAs)) {
    n_na <- sum(is.na(CT.matrix))
    if (n_na > 0) {
      message(sprintf(
        "[NA imputation] %d NA value(s) replaced with values near the minimum. ",
        n_na
      ), "Check raw data before interpretation.")
      set.seed(3)
      mn <- min(CT.matrix, na.rm = TRUE)
      CT.matrix[is.na(CT.matrix)] <- runif(n_na, mn - 1e-8, mn)
    }
  }

  return(CT.matrix)
}

# Backward-compatible alias (deprecated)
extract.nomalised.CTs <- function(...) {
  .Deprecated("extract.normalised.CTs",
              msg = paste("extract.nomalised.CTs() is deprecated.",
                          "Use extract.normalised.CTs()."))
  extract.normalised.CTs(...)
}

# ==============================================================================
# extract.CTs — returns raw CT values (no normalisation)
# ==============================================================================
extract.CTs <- function(CTs, undetermined = 40) {
  CTs <- CTs[, c("Target.Name", "Sample.Name", "CT")]
  CTs[which(CTs$CT == "Undetermined"), "CT"] <- undetermined
  CTs$CT <- gsub(",", ".", CTs$CT)
  CTs$CT <- as.numeric(CTs$CT)

  not.in.use <- c("RT", "RQ1", "PCR", "GDNA", "gDNA", "CTRN", "RQ2")
  CTs <- CTs[which(!(CTs$Target.Name %in% not.in.use)), ]

  samples <- unique(CTs$Sample.Name)

  CT.matrix <- lapply(samples, function(sample) {
    columns <- CTs[which(CTs$Sample.Name == sample), ]
    rownames(columns) <- columns$Target.Name
    columns$Target.Name <- NULL
    colnames(columns)[which(colnames(columns) == "CT")] <- sample
    columns$Sample.Name <- NULL
    return(columns)
  })

  CT.matrix <- do.call("cbind", CT.matrix)
  CT.matrix <- as.matrix(CT.matrix[order(row.names(CT.matrix)), ])
  return(CT.matrix)
}
