# ==============================================================================
# Predefined gene signatures — Stassi Lab
# ==============================================================================

# Metabolism
sign_metabolism <- c('ABCA1', 'ABCG2','ADORA1','AGTR1','AKR1C2','APOA4','CA9',
                     'ENO1','HIF1A', 'LDHB','NAMPT','NR3C1','PDK3','PKM','PRKAA2',
                     'SIRT2','ABCB1', 'ACE','ADRB2','AHR','ALPL','APOA5','CAT',
                     'ESR1','HMOX1','LEP', 'NOS2','NR3C2','PDK4','PPARA','PRKAB1',
                     'SIRT3','ABCC3','ADIPOQ', 'AGER','AIFM1','APOA1','APOD',
                     'CRP','GDF15','HNF1A','LEPR','NOS3', 'PDK1','PIK3CA','PPARG',
                     'PRKACB','TSHR','ABCC5','ADM','AGT', 'AKR1B10','APOA2','AR',
                     'DNMT1','GSTP1','LDHA','MTOR','NOX4','PDK2', 'PIK3CB','PPARGC1',
                     'SIRT1')

# EMT
sign_emt <- c('ACTA2','AEBP1','CLDN1','COL4A1','FSTL1','HMGA2','NOTCH2','SNAI2',
              'TGFBR2','WNT2','ZEB2','ACTB','CALD1','COL1A1','COL4A2','FN1','JAG1',
              'RHOA','TGFB1','TWIST1','WNT5A','ACTG2','CDH1','COL1A2','COL5A2',
              'GREM1','LOX','SMAD4','TGFB2','VIM','WNT7A','ADAM17','CDH2','COL3A1',
              'CTNNB1','HMGA1','NOTCH1','SNAI1','TGFBR1','WNT1','ZEB1')

# Metastasis
sign_metastasis <- c('ANGPT2','APOE','CCN1','CSF1','EGF','ERBB3','FGF2','FGFR4',
                     'IGF1R','IGFBP7','LGALS1','MMP1','MMP9','PDGFRA','PTK2',
                     'SERPINE1','VCAM1','ANGPTL4','APP','CCN2','CSF1R','EGFR',
                     'ETS1','FGFR1','FLT1','IGF2','KDR','LGALS3','MMP2','NCAM1',
                     'PDGFRB','S100A4','SPARC','VEGFA','ANXA1','CALU','CD44',
                     'CTSB','ENG','F2','FGFR2','HGF','IGFBP2','KIT','MET','MMP3',
                     'PDGFA','PLK1','S100A9','SPP1','VEGFB','ANXA2','CAV1','CD68',
                     'EDN1','ERBB2','F2R','FGFR3','IGF1','IGFBP3','LCN2','MIF',
                     'MMP7','PDGFB','POSTN','S100B','TNC','VEGFC')

# List of predefined signatures (used in the analysis pipeline)
lista_firme <- list(
  all                  = "all",
  signature_metabolism = sign_metabolism,
  signature_emt        = sign_emt,
  signature_metastasis = sign_metastasis
)

# ==============================================================================
# get_signatures(cfg)
# ------------------------------------------------------------------------------
# Returns the full list of gene signatures, merging predefined signatures with
# any custom signatures defined in the YAML under the 'custom_signatures' key.
#
# Example YAML:
#   custom_signatures:
#     MySignature: ["GENE1", "GENE2", "GENE3"]
#     AnotherSig:  ["EGFR", "MET", "HGF"]
#
# Custom signatures are stored with key "custom_<name>" and displayed with
# label "Custom: <name>" in the UI.
# ==============================================================================
get_signatures <- function(cfg = NULL) {
  all_sigs <- lista_firme

  if (!is.null(cfg) && !is.null(cfg$custom_signatures) &&
      length(cfg$custom_signatures) > 0) {
    for (sig_name in names(cfg$custom_signatures)) {
      genes <- cfg$custom_signatures[[sig_name]]
      if (is.character(genes) && length(genes) > 0) {
        key <- paste0("custom_", sig_name)
        all_sigs[[key]] <- genes
        message(sprintf("[Custom signature] '%s' loaded: %d genes.", sig_name, length(genes)))
      } else {
        warning(sprintf("[Custom signature] '%s' ignored: must be a non-empty character vector.", sig_name))
      }
    }
  }

  all_sigs
}

# ==============================================================================
# build_signature_choices(cfg)
# ------------------------------------------------------------------------------
# Returns a named vector for selectInput() with all available signatures,
# including any custom ones defined in the YAML.
# ==============================================================================
build_signature_choices <- function(cfg = NULL) {
  base_choices <- c(
    "All genes"  = "all",
    "EMT"        = "signature_emt",
    "Metastasis" = "signature_metastasis",
    "Metabolism" = "signature_metabolism"
  )

  if (!is.null(cfg) && !is.null(cfg$custom_signatures) &&
      length(cfg$custom_signatures) > 0) {
    custom_names   <- names(cfg$custom_signatures)
    custom_vals    <- paste0("custom_", custom_names)
    custom_labels  <- paste0("Custom: ", custom_names)
    custom_choices <- setNames(custom_vals, custom_labels)
    return(c(base_choices, custom_choices))
  }

  base_choices
}
