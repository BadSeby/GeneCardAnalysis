# ==============================================================================
# Tests for extract.normalised.CTs()
# ==============================================================================

# ── Output structure ───────────────────────────────────────────────────────────

test_that("returns a numeric matrix", {
  mat <- extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  expect_true(is.matrix(mat))
  expect_true(is.numeric(mat))
})

test_that("dimensions are genes x samples (housekeeping excluded)", {
  mat <- extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  # 2 target genes (GENE1, GENE2), 2 samples (UT, EMT)
  expect_equal(nrow(mat), 2)
  expect_equal(ncol(mat), 2)
})

test_that("housekeeping gene is absent from row names", {
  mat <- extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  expect_false("GAPDH" %in% rownames(mat))
})

test_that("target genes are present as row names", {
  mat <- extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  expect_true("GENE1" %in% rownames(mat))
  expect_true("GENE2" %in% rownames(mat))
})

test_that("row names are in alphabetical order", {
  mat <- extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  expect_equal(rownames(mat), sort(rownames(mat)))
})

test_that("column names match sample names", {
  mat <- extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  expect_setequal(colnames(mat), c("UT", "EMT"))
})

# ── Mathematical correctness (delta-CT method) ─────────────────────────────────

test_that("2^(-deltaCT) values are mathematically correct", {
  mat <- extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  # UT: GAPDH=20, GENE1=25 → deltaCT=5 → 2^(-5)=0.03125
  expect_equal(mat["GENE1", "UT"], 2^(-5), tolerance = 1e-9)
  # UT: GAPDH=20, GENE2=30 → deltaCT=10 → 2^(-10)
  expect_equal(mat["GENE2", "UT"], 2^(-10), tolerance = 1e-9)
  # EMT: GAPDH=20, GENE1=23 → deltaCT=3 → 2^(-3)=0.125
  expect_equal(mat["GENE1", "EMT"], 2^(-3), tolerance = 1e-9)
  # EMT: GAPDH=20, GENE2=28 → deltaCT=8 → 2^(-8)
  expect_equal(mat["GENE2", "EMT"], 2^(-8), tolerance = 1e-9)
})

test_that("all output values are strictly positive", {
  mat <- extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  expect_true(all(mat > 0, na.rm = TRUE))
})

# ── "Undetermined" handling ────────────────────────────────────────────────────

test_that("'Undetermined' string is replaced with undetermined value", {
  df <- make_cts_undetermined()
  # With handle.NAs = FALSE so we can inspect the raw replacement
  # GENE2/UT becomes Undetermined → replaced by 40 → deltaCT = 40 - 20 = 20 → 2^(-20)
  mat <- suppressMessages(
    extract.normalised.CTs(df, housekeeping.gene = "GAPDH",
                           handle.NAs = FALSE, undetermined = 40)
  )
  expect_equal(mat["GENE2", "UT"], 2^(-20), tolerance = 1e-9)
})

test_that("custom undetermined value is respected", {
  df <- make_cts_undetermined()
  mat <- suppressMessages(
    extract.normalised.CTs(df, housekeeping.gene = "GAPDH",
                           handle.NAs = FALSE, undetermined = 35)
  )
  # GENE2/UT: CT=35, GAPDH=20 → deltaCT=15 → 2^(-15)
  expect_equal(mat["GENE2", "UT"], 2^(-15), tolerance = 1e-9)
})

# ── Comma decimal separator ────────────────────────────────────────────────────

test_that("CT values with comma as decimal separator are parsed correctly", {
  mat_comma  <- suppressMessages(
    extract.normalised.CTs(make_cts_comma(), housekeeping.gene = "GAPDH")
  )
  mat_period <- suppressMessages(
    extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  )
  expect_equal(mat_comma, mat_period, tolerance = 1e-9)
})

# ── Exclude genes ──────────────────────────────────────────────────────────────

test_that("control genes are excluded from the output", {
  mat <- suppressMessages(
    extract.normalised.CTs(make_cts_with_control(), housekeeping.gene = "GAPDH")
  )
  expect_false("gDNA" %in% rownames(mat))
})

test_that("housekeeping gene is NOT removed by exclude_genes", {
  # GAPDH is in both exclude_genes and housekeeping.gene — must survive for normalisation
  mat <- suppressMessages(
    extract.normalised.CTs(
      make_cts_simple(),
      housekeeping.gene = "GAPDH",
      exclude_genes     = c("GAPDH", "RT")   # GAPDH explicitly listed
    )
  )
  # Must still be absent from output (removed after use), but math must be correct
  expect_false("GAPDH" %in% rownames(mat))
  expect_equal(mat["GENE1", "UT"], 2^(-5), tolerance = 1e-9)
})

# ── Duplicate averaging ────────────────────────────────────────────────────────

test_that("technical duplicates are averaged before normalisation", {
  mat_dup    <- suppressMessages(
    extract.normalised.CTs(make_cts_duplicates(), housekeeping.gene = "GAPDH")
  )
  mat_simple <- suppressMessages(
    extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  )
  # Duplicates are (orig + orig+delta)/2; result should differ slightly from simple
  # but remain close (delta = 0.1–0.3 CT units → small effect on 2^(-deltaCT))
  expect_true(all(abs(mat_dup - mat_simple) < 0.05))
  # Dimensions unchanged
  expect_equal(dim(mat_dup), dim(mat_simple))
})

test_that("partial duplicates trigger a warning", {
  expect_warning(
    suppressMessages(
      extract.normalised.CTs(make_cts_partial_dup(), housekeeping.gene = "GAPDH")
    ),
    regexp = "Not all genes are duplicated"
  )
})

# ── Multiple housekeeping genes ────────────────────────────────────────────────

test_that("multiple housekeeping genes are averaged for delta-CT", {
  mat <- suppressMessages(
    extract.normalised.CTs(make_cts_two_hk(), housekeeping.gene = c("GAPDH", "HPRT1"))
  )
  # UT: mean(GAPDH=20, HPRT1=21) = 20.5
  # GENE1/UT: deltaCT = 25 - 20.5 = 4.5 → 2^(-4.5)
  expect_equal(mat["GENE1", "UT"], 2^(-4.5), tolerance = 1e-9)
  # Both HK genes absent from output
  expect_false("GAPDH" %in% rownames(mat))
  expect_false("HPRT1" %in% rownames(mat))
})

test_that("warning is issued when one of two HK genes is missing from a sample", {
  df <- make_cts_simple()   # only GAPDH, no HPRT1
  expect_warning(
    suppressMessages(
      extract.normalised.CTs(df, housekeeping.gene = c("GAPDH", "HPRT1"))
    ),
    regexp = "Missing housekeeping genes"
  )
})

# ── Error on absent housekeeping gene ─────────────────────────────────────────

test_that("stops with error when no housekeeping gene is found in a sample", {
  expect_error(
    suppressMessages(
      extract.normalised.CTs(make_cts_no_hk(), housekeeping.gene = "GAPDH")
    ),
    regexp = "No housekeeping gene found"
  )
})

# ── Auto-selection of housekeeping gene ───────────────────────────────────────

test_that("auto-selects HPRT1 when it has lower variance than GAPDH", {
  # make_cts_autoselect_hprt1: HPRT1 is constant at 21 across all samples
  expect_message(
    extract.normalised.CTs(make_cts_autoselect_hprt1()),
    regexp = "HPRT1.*selected automatically"
  )
})

test_that("auto-selects GAPDH when it has lower variance than HPRT1", {
  df <- data.frame(
    Sample.Name = rep(c("S1", "S2"), each = 3),
    Target.Name = rep(c("GAPDH", "HPRT1", "GENE1"), times = 2),
    CT          = c(20, 21, 25,   20, 24, 26),  # GAPDH var=0, HPRT1 var=4.5
    stringsAsFactors = FALSE
  )
  expect_message(
    extract.normalised.CTs(df),
    regexp = "GAPDH.*selected automatically"
  )
})

# ── NA imputation ──────────────────────────────────────────────────────────────

test_that("NAs are imputed when handle.NAs = TRUE", {
  df <- make_cts_simple()
  df$CT[df$Sample.Name == "UT" & df$Target.Name == "GENE2"] <- NA
  df$CT <- as.character(df$CT)
  mat <- suppressMessages(
    extract.normalised.CTs(df, housekeeping.gene = "GAPDH", handle.NAs = TRUE)
  )
  expect_false(any(is.na(mat)))
})

test_that("NAs are preserved when handle.NAs = FALSE", {
  df <- make_cts_simple()
  df$CT[df$Sample.Name == "UT" & df$Target.Name == "GENE2"] <- NA
  df$CT <- as.character(df$CT)
  mat <- suppressMessages(
    extract.normalised.CTs(df, housekeeping.gene = "GAPDH", handle.NAs = FALSE)
  )
  expect_true(any(is.na(mat)))
})

test_that("NA imputation emits an informative message", {
  df <- make_cts_simple()
  df$CT[1] <- NA
  df$CT <- as.character(df$CT)
  expect_message(
    extract.normalised.CTs(df, housekeeping.gene = "GAPDH", handle.NAs = TRUE),
    regexp = "NA imputation"
  )
})

# ── Deprecated alias ──────────────────────────────────────────────────────────

test_that("extract.nomalised.CTs() is a working deprecated alias", {
  mat_new  <- suppressMessages(
    extract.normalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  )
  mat_old  <- suppressWarnings(suppressMessages(
    extract.nomalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
  ))
  expect_equal(mat_new, mat_old)
})

test_that("extract.nomalised.CTs() emits a deprecation warning", {
  expect_warning(
    suppressMessages(
      extract.nomalised.CTs(make_cts_simple(), housekeeping.gene = "GAPDH")
    ),
    regexp = "deprecated"
  )
})

# ==============================================================================
# Tests for extract.CTs()
# ==============================================================================

test_that("extract.CTs returns a numeric matrix", {
  mat <- extract.CTs(make_cts_simple())
  expect_true(is.matrix(mat))
  expect_true(is.numeric(mat))
})

test_that("extract.CTs returns raw CT values (not 2^-deltaCT)", {
  mat <- extract.CTs(make_cts_simple())
  # GAPDH/UT should be 20, not 2^(-something)
  expect_equal(mat["GAPDH", "UT"], 20)
})

test_that("extract.CTs replaces 'Undetermined' with the numeric value", {
  df <- make_cts_undetermined()
  mat <- extract.CTs(df, undetermined = 40)
  expect_equal(mat["GENE2", "UT"], 40)
})

test_that("extract.CTs excludes standard control genes", {
  df <- make_cts_with_control()
  mat <- extract.CTs(df)
  expect_false("gDNA" %in% rownames(mat))
})

test_that("extract.CTs handles comma decimal separators", {
  mat_comma  <- extract.CTs(make_cts_comma())
  mat_period <- extract.CTs(make_cts_simple())
  expect_equal(mat_comma, mat_period, tolerance = 1e-9)
})

test_that("extract.CTs row names are in alphabetical order", {
  mat <- extract.CTs(make_cts_simple())
  expect_equal(rownames(mat), sort(rownames(mat)))
})
