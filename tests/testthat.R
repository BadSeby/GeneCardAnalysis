library(testthat)

# Source the package functions directly (pre-package structure)
source(file.path(dirname(dirname(getwd())), "card_functions.R"), local = TRUE)
source(file.path(dirname(dirname(getwd())), "signatures_genes_card.R"), local = TRUE)

test_check("GeneCardAnalysis")
