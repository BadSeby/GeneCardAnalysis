# ==============================================================================
# Tests for get_signatures() and build_signature_choices()
# ==============================================================================

# ── get_signatures() ──────────────────────────────────────────────────────────

test_that("get_signatures() without cfg returns the four base signatures", {
  sigs <- get_signatures()
  expect_true(is.list(sigs))
  expect_true("all"                  %in% names(sigs))
  expect_true("signature_emt"        %in% names(sigs))
  expect_true("signature_metastasis" %in% names(sigs))
  expect_true("signature_metabolism" %in% names(sigs))
})

test_that("'all' signature is the string 'all'", {
  sigs <- get_signatures()
  expect_equal(sigs$all, "all")
})

test_that("built-in signatures are non-empty character vectors", {
  sigs <- get_signatures()
  for (nm in c("signature_emt", "signature_metastasis", "signature_metabolism")) {
    expect_true(is.character(sigs[[nm]]))
    expect_gt(length(sigs[[nm]]), 0)
  }
})

test_that("custom signature is added with 'custom_' prefix", {
  cfg <- list(custom_signatures = list(MyPanel = c("EGFR", "MET", "HGF")))
  sigs <- suppressMessages(get_signatures(cfg))
  expect_true("custom_MyPanel" %in% names(sigs))
  expect_equal(sigs$custom_MyPanel, c("EGFR", "MET", "HGF"))
})

test_that("multiple custom signatures are all added", {
  cfg <- list(custom_signatures = list(
    PanelA = c("GENE1", "GENE2"),
    PanelB = c("GENE3", "GENE4", "GENE5")
  ))
  sigs <- suppressMessages(get_signatures(cfg))
  expect_true("custom_PanelA" %in% names(sigs))
  expect_true("custom_PanelB" %in% names(sigs))
  expect_equal(length(sigs$custom_PanelB), 3)
})

test_that("custom signature loading emits an informative message", {
  cfg <- list(custom_signatures = list(TestSig = c("A", "B")))
  expect_message(get_signatures(cfg), regexp = "TestSig.*loaded.*2 genes")
})

test_that("invalid custom signature (empty vector) is ignored with a warning", {
  cfg <- list(custom_signatures = list(BadSig = character(0)))
  expect_warning(get_signatures(cfg), regexp = "BadSig.*ignored")
  sigs <- suppressWarnings(get_signatures(cfg))
  expect_false("custom_BadSig" %in% names(sigs))
})

test_that("invalid custom signature (non-character) is ignored with a warning", {
  cfg <- list(custom_signatures = list(NumSig = c(1, 2, 3)))
  expect_warning(get_signatures(cfg), regexp = "NumSig.*ignored")
  sigs <- suppressWarnings(get_signatures(cfg))
  expect_false("custom_NumSig" %in% names(sigs))
})

test_that("NULL cfg returns only base signatures", {
  expect_equal(get_signatures(NULL), get_signatures())
})

test_that("cfg with no custom_signatures key returns only base signatures", {
  cfg <- list(project_name = "test", comparisons = list())
  expect_equal(get_signatures(cfg), get_signatures())
})

test_that("base signatures are unchanged after adding custom ones", {
  cfg <- list(custom_signatures = list(Extra = c("X1", "X2")))
  sigs_base   <- get_signatures()
  sigs_custom <- suppressMessages(get_signatures(cfg))
  for (nm in names(sigs_base)) {
    expect_equal(sigs_custom[[nm]], sigs_base[[nm]])
  }
})

# ── build_signature_choices() ─────────────────────────────────────────────────

test_that("returns a named character vector", {
  choices <- build_signature_choices()
  expect_true(is.character(choices))
  expect_false(is.null(names(choices)))
})

test_that("returns exactly 4 base choices without cfg", {
  choices <- build_signature_choices()
  expect_equal(length(choices), 4)
})

test_that("base choice values are correct", {
  choices <- build_signature_choices()
  expect_true("all"                  %in% choices)
  expect_true("signature_emt"        %in% choices)
  expect_true("signature_metastasis" %in% choices)
  expect_true("signature_metabolism" %in% choices)
})

test_that("base choice labels are human-readable", {
  choices <- build_signature_choices()
  expect_true("All genes"  %in% names(choices))
  expect_true("EMT"        %in% names(choices))
  expect_true("Metastasis" %in% names(choices))
  expect_true("Metabolism" %in% names(choices))
})

test_that("custom signatures appear as additional choices", {
  cfg <- list(custom_signatures = list(MyPanel = c("A", "B")))
  choices <- build_signature_choices(cfg)
  expect_equal(length(choices), 5)
  expect_true("custom_MyPanel" %in% choices)
  expect_true("Custom: MyPanel" %in% names(choices))
})

test_that("multiple custom signatures each get a labelled choice entry", {
  cfg <- list(custom_signatures = list(P1 = c("A"), P2 = c("B", "C")))
  choices <- build_signature_choices(cfg)
  expect_equal(length(choices), 6)
  expect_true("Custom: P1" %in% names(choices))
  expect_true("Custom: P2" %in% names(choices))
})

test_that("NULL cfg returns the same 4 base choices as no-arg call", {
  expect_equal(build_signature_choices(NULL), build_signature_choices())
})

test_that("cfg without custom_signatures returns only base choices", {
  cfg <- list(project_name = "test")
  expect_equal(build_signature_choices(cfg), build_signature_choices())
})
