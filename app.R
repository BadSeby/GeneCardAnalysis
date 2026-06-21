# ==============================================================================
# Gene Card Analysis — Shiny App
# Stassi Lab
# RT-qPCR data analysis: delta-CT normalisation, heatmaps, log2FC
# ==============================================================================

library(shiny)
library(bslib)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(tidyverse)
library(grid)
library(glue)
library(writexl)
library(showtext)
library(yaml)
library(DT)
library(colourpicker)
library(gtools)

source("card_functions.R")
source("signatures_genes_card.R")

# Font and showtext — must be called once at startup
tryCatch(
  font_add_google("Arimo", "Arimo"),
  error = function(e) {
    message("Google Font unavailable (", conditionMessage(e), "). Using default font.")
  }
)
showtext_auto()

# ==============================================================================
# Base signature choices (fallback; dynamic choices are in reactive_signature_choices in server)
signature_choices_base <- c(
  "All genes"  = "all",
  "EMT"        = "signature_emt",
  "Metastasis" = "signature_metastasis",
  "Metabolism" = "signature_metabolism"
)

# Null-coalescing helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# HELPER: plot height in pixels (used in both renderPlot and download handlers)
# ==============================================================================
plot_h_px <- function(n_genes, sig_type = "all", dpi = 96) {
  cell_h_mm <- if (sig_type == "all") max(3, min(5, 400 / n_genes))
  else                   max(8, 300 / n_genes)
  h_in <- (n_genes * cell_h_mm / 25.4) + 4
  max(400L, round(h_in * dpi))
}

# ==============================================================================
# HELPER: draw ComplexHeatmap on current device
# font_mult = 1 (screen/PDF), font_mult = 3 (high-resolution PNG)
# ==============================================================================
draw_heatmap_cht <- function(mat, col_fun, legend_name, legend_title,
                             comp_display_name, comp_info, cfg, author_name,
                             font_mult = 1, interp_override = NULL,
                             cluster_override = NULL) {
  n_genes  <- nrow(mat)
  sig_type <- comp_info$signature_type %||% "all"
  # Determine row/column clustering — UI override takes priority over YAML
  cluster_type <- if (!is.null(cluster_override)) {
    cluster_override
  } else {
    comp_info$cluster_type %||%
      if (isTRUE(comp_info$cluster_heatmap)) "both" else "none"
  }
  do_clust_r <- cluster_type %in% c("rows", "both")
  do_clust_c <- cluster_type %in% c("cols", "both")
  ut_ref   <- comp_info$reference %||% "—"
  
  # Adaptive cell height
  cell_h_mm    <- if (sig_type == "all") max(3, min(5, 400 / n_genes))
  else                   max(8, 300 / n_genes)
  row_fs       <- max(6, min(12, cell_h_mm * 1.5))
  heatmap_h_mm <- n_genes * cell_h_mm

  # Annotation text
  gs      <- cfg$global_settings %||% list()
  hk_used <- paste(gs$housekeeping_genes %||% "—", collapse = ", ")

  meta_txt <- paste0(
    "Summary Statistics:\n",
    "- Genes: ",          n_genes,      "\n",
    "- Samples: ",        ncol(mat),    "\n",
    "- Housekeeping: ",   hk_used,      "\n",
    "- Reference (UT): ", ut_ref
  )
  interp_txt <- if (!is.null(interp_override)) {
    interp_override
  } else if (grepl("z.?score|zscore|z_score", legend_name, ignore.case = TRUE)) {
    paste0("Z-score interpretation:\n",
           "> 0 = above gene mean\n",
           "< 0 = below gene mean\n",
           "= 0 = at gene mean\n",
           "Unit = standard deviations (log2 scale)")
  } else {
    paste0("Log2FC interpretation:\n",
           "> 0 = upregulated\n",
           "< 0 = downregulated\n",
           "= 0 = unchanged\n",
           "Each unit = 2-fold difference")
  }
  final_meta  <- paste0(meta_txt, "\n\n", interp_txt)
  caption_txt <- paste0(
    "Stassi Lab | Project: ", cfg$project_name %||% "—",
    "\nAuthor: ", author_name, " | Date: ", Sys.Date()
  )
  
  ht <- Heatmap(
    mat,
    name             = legend_name,
    width            = unit(50, "mm"),
    height           = unit(heatmap_h_mm, "mm"),
    col              = col_fun,
    cluster_columns  = do_clust_c,
    cluster_rows     = do_clust_r,
    rect_gp          = gpar(col = "white", lwd = 0.6 * font_mult),
    column_title     = comp_display_name,
    column_title_gp  = gpar(fontsize = 16 * font_mult, fontface = "bold"),
    row_names_gp     = gpar(fontsize = row_fs * font_mult),
    column_names_gp  = gpar(fontsize = 11 * font_mult, fontface = "bold"),
    column_names_rot = 45,
    heatmap_legend_param = list(
      title     = legend_title,
      title_gp  = gpar(fontsize = 10 * font_mult, fontface = "bold"),
      labels_gp = gpar(fontsize =  8 * font_mult)
    )
  )
  
  draw(ht,
       heatmap_legend_side = "right",
       padding = unit(c(25, 25, 25, 25), "mm"))
  
  # Fixed font size for annotation text (not scaled with font_mult)
  txt_fs_cap <- if (font_mult > 1) 24 else 8
  txt_fs_met <- if (font_mult > 1) 21 else 7
  lh         <- if (font_mult > 1) 0.5 else 1.1
  
  grid.text(caption_txt,
            x = unit(0.98, "npc"), y = unit(0.01, "npc"),
            just = c("right", "bottom"),
            gp   = gpar(fontsize = txt_fs_cap, col = "grey30", lineheight = lh))
  
  grid.text(final_meta,
            x = unit(0.02, "npc"), y = unit(0.01, "npc"),
            just = c("left", "bottom"),
            gp   = gpar(fontsize = txt_fs_met, col = "black",  lineheight = lh))
}

# ==============================================================================
# HELPER: ranked bar chart for Log2FC (used when <= 2 samples are present)
# Shows all genes ordered by log2FC, red bars (up) / blue bars (down).
# ==============================================================================
draw_log2fc_bars <- function(log2fc_df, comp_info, cfg, author_name) {
  ref_col     <- comp_info$reference %||% ""
  sample_cols <- setdiff(colnames(log2fc_df)[-1], ref_col)
  if (length(sample_cols) == 0) sample_cols <- colnames(log2fc_df)[-1]
  samp <- sample_cols[1]

  df <- data.frame(
    Gene   = log2fc_df$Gene,
    log2FC = as.numeric(log2fc_df[[samp]]),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$log2FC), ]
  df <- df[order(df$log2FC), ]
  df$Gene      <- factor(df$Gene, levels = df$Gene)
  df$Direction <- ifelse(df$log2FC >= 0, "Up", "Down")

  caption_txt <- paste0(
    "Stassi Lab | Project: ", cfg$project_name %||% "—",
    "\nAuthor: ", author_name, " | Date: ", Sys.Date()
  )

  ggplot(df, aes(x = log2FC, y = Gene, fill = Direction)) +
    geom_col(width = 0.72, alpha = 0.88) +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey20") +
    scale_fill_manual(
      values = c("Up" = "#b2182b", "Down" = "#2166ac"),
      guide  = guide_legend(title = NULL)
    ) +
    scale_x_continuous(expand = expansion(mult = 0.06)) +
    labs(
      title    = paste0(comp_info$name %||% "", "  ·  Log2FC vs. ", ref_col),
      subtitle = paste0("Sample: ", samp, "  ·  Genes: ", nrow(df)),
      caption  = caption_txt,
      x        = "Log2 Fold Change",
      y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = rel(0.72)),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(colour = "grey40", size = 10),
      plot.caption       = element_text(colour = "grey50", size = 7, hjust = 0),
      legend.position    = "top"
    )
}

# ==============================================================================
# UI — Gene Card Analysis (dark sidebar, icone Font Awesome)
# Richiede il file  www/custom.css  nella cartella dell'app.
# ==============================================================================
ui <- page_sidebar(
  title = tags$span(
    class = "app-brand",
    tags$span(class = "app-logo", icon("dna")),
    tags$span(
      class = "app-brand-text",
      tags$span(class = "app-brand-title", "Gene Card Analysis"),
      tags$span(class = "app-brand-sub", "Stassi Lab")
    )
  ),
  window_title = "Gene Card Analysis — Stassi Lab",
  
  theme = bs_theme(
    version    = 5,
    bootswatch = "flatly",
    primary    = "#3b82f6",
    secondary  = "#8b5cf6",
    success    = "#10b981",
    warning    = "#f59e0b",
    danger     = "#ef4444",
    info       = "#06b6d4",
    base_font  = font_google("Arimo"),
    "border-radius" = "0.7rem"
  ),
  
  tags$head(
    tags$link(rel = "stylesheet", href = "custom.css?v=4")
  ),
  
  # ── SIDEBAR ────────────────────────────────────────────────────────────────
  sidebar = sidebar(
    width = 330,
    class = "gca-sidebar",
    
    div(class = "side-label", icon("folder-open"), "Load Data"),
    card(
      class = "gca-side-card",
      fileInput("file_csv", "CT file (.csv / .txt, tab-separated)",
                accept = c(".csv", ".txt", ".tsv"),
                buttonLabel = "Browse…",
                placeholder = "No file selected"),
      fileInput("file_yaml", "Configuration (.yaml / .txt)",
                accept = c(".yaml", ".yml", ".txt",
                           "text/yaml", "text/x-yaml",
                           "application/x-yaml", "text/plain"),
                buttonLabel = "Browse…",
                placeholder = "No file selected"),
      radioButtons("config_mode", label = "Configuration source:",
                   choices  = c("From loaded YAML" = "yaml",
                                "Manual form"       = "form"),
                   selected = "yaml",
                   inline   = TRUE)
    ),

    uiOutput("yaml_load_status"),

    div(class = "side-label", icon("sliders"), "Global Settings"),
    card(
      class = "gca-side-card",
      textInput("project_name", "Project Name", value = ""),
      textInput("author_name", "Author", value = "",
                placeholder = "Your name"),
      numericInput("undetermined_value", "Undetermined Value", value = 40),
      textInput("housekeeping_genes", "Housekeeping Genes (comma-separated)",
                value = "GAPDH, HPRT1"),
      textInput("exclude_genes", "Genes to exclude (comma-separated)",
                value = "RT, gDNA")
    ),

    div(class = "side-label", icon("palette"), "Heatmap Appearance"),
    card(
      class = "gca-side-card",
      div(
        style = "display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 6px;",
        div(colourInput("col_low",  "Min", "#4575b4", returnName = TRUE)),
        div(colourInput("col_mid",  "Mid", "#f7f7f7", returnName = TRUE)),
        div(colourInput("col_high", "Max", "#d73027", returnName = TRUE))
      ),
      helpText("Applied to all heatmaps (z-score, normalised, log2FC).")
    ),

    uiOutput("form_comparisons_panel"),

    actionButton("btn_run",
                 label = tagList(icon("play"), "Run Analysis"),
                 class = "btn-primary w-100 mt-2 gca-btn-run"),
    uiOutput("yaml_export_btn"),

    hr(),
    helpText("Select 'Manual form' to configure comparisons without a YAML file.")
  ),
  
  # ── MAIN PANEL ─────────────────────────────────────────────────────────────
  navset_card_underline(
    id = "main_tabs",
    
    nav_panel(
      "Home",
      icon = icon("house"),
      layout_columns(
        col_widths = c(8, 4),
        card(
          card_header(icon("circle-info"), "Welcome to Gene Card Analysis"),
          p("This app analyses RT-qPCR CT values and generates comparative heatmaps."),
          p("Getting started:"),
          tags$ol(
            tags$li("Load the CT data file (tab-separated CSV with columns: Sample.Name, Target.Name, CT)"),
            tags$li("Load a YAML configuration file or fill in the manual form in the sidebar"),
            tags$li("Click ", tags$b("Run Analysis"))
          ),
          hr(),
          p(tags$b("Outputs generated per comparison:")),
          tags$ul(
            tags$li("Z-score heatmap (interactive view + PNG/PDF download)"),
            tags$li("Normalised values heatmap (PDF download)"),
            tags$li("Interactive log2FC table (Excel download)"),
            tags$li("Log2FC heatmap / ranked bar chart (PDF download)")
          )
        ),
        card(
          card_header(icon("signal"), "Status"),
          uiOutput("status_box")
        )
      )
    ),
    
    nav_panel(
      "Data",
      icon = icon("table"),
      card(
        card_header(icon("table-list"), "Loaded CT file preview"),
        uiOutput("data_summary"),
        DTOutput("table_preview")
      )
    ),

    nav_panel(
      "Configuration",
      icon = icon("sliders"),
      card(
        card_header(icon("layer-group"), "Defined comparisons"),
        uiOutput("config_summary"),
        uiOutput("comparisons_ui")
      )
    ),

    nav_panel(
      "Results",
      icon = icon("chart-simple"),
      uiOutput("results_ui")
    ),

    nav_panel(
      "Download",
      icon = icon("download"),
      uiOutput("downloads_ui")
    ),

    nav_panel(
      "Log",
      icon = icon("terminal"),
      card(
        card_header(icon("terminal"), "Execution log"),
        verbatimTextOutput("log_output")
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================
server <- function(input, output, session) {
  
  # ── Application state ───────────────────────────────────────────────────────
  rv <- reactiveValues(
    log_lines    = character(0),
    results      = list(),
    comp_ids     = integer(0),
    comp_counter = 0L,
    comp_presets = list(),
    yaml_error   = NULL
  )
  
  log_msg <- function(msg, type = "info") {
    prefix <- switch(type, info = "ℹ", ok = "✅", warn = "⚠️", error = "❌", "ℹ")
    ts <- format(Sys.time(), "[%H:%M:%S]")
    rv$log_lines <- c(isolate(rv$log_lines), paste(ts, prefix, msg))
  }
  
  # ── REACTIVE: CSV reading ───────────────────────────────────────────────────
  reactive_data <- reactive({
    req(input$file_csv)
    path <- input$file_csv$datapath
    
    df <- tryCatch({
      df_try <- read.csv(path, stringsAsFactors = FALSE, sep = "\t")
      if (ncol(df_try) >= 3) {
        df_try
      } else {
        df_comma <- read.csv(path, stringsAsFactors = FALSE, sep = ",")
        if (ncol(df_comma) >= 3) {
          log_msg("CSV separator auto-detected: comma (,)", "info")
          df_comma
        } else {
          df_semi <- read.csv(path, stringsAsFactors = FALSE, sep = ";")
          if (ncol(df_semi) >= 3) {
            log_msg("CSV separator auto-detected: semicolon (;)", "info")
            df_semi
          } else {
            log_msg("Cannot detect file separator. Expected at least 3 fields (Sample.Name, Target.Name, CT).", "error")
            NULL
          }
        }
      }
    }, error = function(e) {
      log_msg(paste("CSV read error:", e$message), "error")
      NULL
    })
    if (is.null(df)) return(NULL)

    required_cols <- c("Sample.Name", "Target.Name", "CT")
    missing_cols  <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
      log_msg(paste("Missing columns in CSV:", paste(missing_cols, collapse = ", ")), "error")
      return(NULL)
    }

    log_msg(paste0("CSV loaded: ", nrow(df), " rows, ",
                   length(unique(df$Sample.Name)), " samples, ",
                   length(unique(df$Target.Name)), " genes"), "ok")
    df
  })
  
  # ── REACTIVE: YAML parsing ──────────────────────────────────────────────────
  reactive_config <- reactive({
    req(input$file_yaml)
    path <- input$file_yaml$datapath
    
    cfg <- tryCatch(
      yaml::read_yaml(path),
      error = function(e) NULL
    )
    
    if (is.null(cfg)) return(NULL)
    if (is.null(cfg$comparisons) || length(cfg$comparisons) == 0) return(NULL)
    
    cfg
  })
  
  # ── OBSERVER: update rv$yaml_error as a separate side effect ────────────────
  observeEvent(input$file_yaml, {
    path <- input$file_yaml$datapath

    cfg <- tryCatch(
      yaml::read_yaml(path),
      error = function(e) {
        rv$yaml_error <- paste("YAML parsing error:", e$message)
        log_msg(rv$yaml_error, "error")
        NULL
      }
    )

    if (is.null(cfg)) return()

    if (is.null(cfg$comparisons) || length(cfg$comparisons) == 0) {
      rv$yaml_error <- "YAML loaded but no comparisons found under the 'comparisons' key"
      log_msg(rv$yaml_error, "error")
      return()
    }

    rv$yaml_error <- NULL
    log_msg(paste0("YAML loaded: project '", cfg$project_name,
                   "', ", length(cfg$comparisons), " comparison(s)"), "ok")
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  # ── YAML inline status (sidebar) ────────────────────────────────────────────
  output$yaml_load_status <- renderUI({
    if (is.null(input$file_yaml)) return(NULL)

    err <- rv$yaml_error
    cfg <- reactive_config()

    if (!is.null(err)) {
      div(class = "alert alert-danger p-1 mt-1",
          style = "font-size:.8em;",
          icon("circle-xmark"), " ", err)
    } else if (!is.null(cfg)) {
      div(class = "alert alert-success p-1 mt-1",
          style = "font-size:.8em;",
          icon("circle-check"), " ",
          strong(cfg$project_name %||% "—"), " — ",
          length(cfg$comparisons), " comparison(s) loaded")
    } else {
      div(class = "alert alert-warning p-1 mt-1",
          style = "font-size:.8em;",
          icon("spinner"), " Loading…")
    }
  })
  
  # ── REACTIVE: resolved author name ──────────────────────────────────────────
  reactive_author <- reactive({
    name <- trimws(input$author_name %||% "")
    if (nchar(name) > 0) return(name)
    cfg <- reactive_config()
    if (!is.null(cfg) && !is.null(cfg$author)) as.character(cfg$author) else "Unknown Author"
  })
  
  # ── REACTIVE: signature choices (base + custom from YAML) ───────────────────
  reactive_signature_choices <- reactive({
    build_signature_choices(reactive_config())
  })
  
  # ── REACTIVE: comparisons from manual form ──────────────────────────────────
  reactive_form_comparisons <- reactive({
    ids <- rv$comp_ids
    if (length(ids) == 0) return(list())
    
    df           <- reactive_data()
    avail_samp   <- if (!is.null(df)) sort(unique(df$Sample.Name)) else character(0)
    csv_loaded   <- length(avail_samp) > 0
    
    lapply(ids, function(id) {
      key      <- as.character(id)
      preset   <- rv$comp_presets[[key]] %||% list()
      
      name_v   <- input[[paste0("comp_name_",    id)]] %||% preset$name    %||% paste0("Comp", id)
      ref_v    <- input[[paste0("comp_ref_",     id)]] %||% preset$reference %||% ""
      sig_v    <- input[[paste0("comp_sig_",     id)]] %||% preset$signature_type %||% "all"
      clust_v  <- isTRUE(input[[paste0("comp_cluster_", id)]])
      
      samp_v <- if (csv_loaded) {
        sel <- input[[paste0("comp_samples_", id)]]
        if (is.null(sel) || length(sel) == 0) "all" else sel
      } else {
        raw  <- input[[paste0("comp_samples_txt_", id)]] %||% ""
        parsed <- trimws(strsplit(raw, ",")[[1]])
        if (length(parsed) == 0 || all(parsed == "")) "all" else parsed
      }
      
      list(
        name            = name_v,
        reference       = ref_v,
        samples         = samp_v,
        signature_type  = sig_v,
        cluster_heatmap = clust_v
      )
    })
  })
  
  # ── REACTIVE: active config (YAML or form, with global settings override) ───
  reactive_active_config <- reactive({
    hk <- trimws(strsplit(input$housekeeping_genes %||% "GAPDH,HPRT1", ",")[[1]])
    ex <- trimws(strsplit(input$exclude_genes      %||% "RT,gDNA",     ",")[[1]])
    gs_form <- list(
      housekeeping_genes = hk,
      undetermined_value = input$undetermined_value %||% 40,
      exclude_genes      = ex
    )
    author_v <- trimws(input$author_name %||% "Unknown")
    
    if (!is.null(input$config_mode) && input$config_mode == "form") {
      comps <- reactive_form_comparisons()
      if (length(comps) == 0) return(NULL)
      proj_name <- if (nchar(trimws(input$project_name %||% "")) > 0)
        input$project_name else "GeneCardAnalysis"
      return(list(
        project_name    = proj_name,
        author          = author_v,
        global_settings = gs_form,
        comparisons     = comps
      ))
    }
    
    cfg <- reactive_config()
    if (is.null(cfg)) return(NULL)
    
    proj <- trimws(input$project_name %||% "")
    if (nchar(proj) > 0) cfg$project_name <- proj
    
    if (length(hk) > 0 && !all(hk == ""))
      cfg$global_settings$housekeeping_genes <- hk
    if (length(ex) > 0 && !all(ex == ""))
      cfg$global_settings$exclude_genes <- ex
    cfg$global_settings$undetermined_value <- input$undetermined_value %||% 40
    
    cfg
  })
  
  # ── FORM COMPARISONS PANEL ──────────────────────────────────────────────────
  # Rendered only when config_mode == "form"
  output$form_comparisons_panel <- renderUI({
    if (is.null(input$config_mode) || input$config_mode != "form") return(NULL)
    
    df         <- reactive_data()
    avail_samp <- if (!is.null(df)) sort(unique(df$Sample.Name)) else character(0)
    csv_loaded <- length(avail_samp) > 0
    ids        <- rv$comp_ids
    
    comp_cards <- lapply(ids, function(id) {
      key    <- as.character(id)
      preset <- rv$comp_presets[[key]] %||% list()
      
      preset_samp <- preset$samples %||% character(0)
      if (identical(preset_samp, "all")) preset_samp <- character(0)
      
      card(
        style = "border-left: 3px solid #3b82f6; margin-bottom:8px;",
        card_header(
          class = "d-flex justify-content-between align-items-center py-1",
          tags$span(style = "font-size:.85em; font-weight:bold;",
                    icon("microscope"), paste0(" Comparison #", id)),
          actionButton(paste0("btn_rm_", id), icon("xmark"),
                       class = "btn-sm btn-outline-danger py-0 px-2")
        ),
        div(style = "padding: 6px 8px;",
            textInput(paste0("comp_name_", id), "Name",
                      value = preset$name %||% paste0("Comp", id),
                      placeholder = "e.g. EMT_vs_UT"),
            if (csv_loaded) {
              tagList(
                selectInput(paste0("comp_ref_", id), "Reference (UT)",
                            choices  = c("— select —" = "", avail_samp),
                            selected = preset$reference %||% ""),
                selectInput(paste0("comp_samples_", id),
                            "Samples (empty = all)",
                            choices  = avail_samp,
                            selected = preset_samp,
                            multiple = TRUE)
              )
            } else {
              tagList(
                textInput(paste0("comp_ref_",     id), "Reference (UT)",
                          value = preset$reference %||% "",
                          placeholder = "exact sample name"),
                textInput(paste0("comp_samples_txt_", id),
                          "Samples (comma-separated, empty = all)",
                          value = if (length(preset_samp) > 0)
                            paste(preset_samp, collapse = ", ")
                          else "",
                          placeholder = "all")
              )
            },
            selectInput(paste0("comp_sig_", id), "Gene signature",
                        choices  = reactive_signature_choices(),
                        selected = preset$signature_type %||% "all"),
            checkboxInput(paste0("comp_cluster_", id),
                          "Cluster rows/columns",
                          value = isTRUE(preset$cluster_heatmap))
        )
      )
    })
    
    btns <- div(
      actionButton("btn_add_comp",
                   label = tagList(icon("plus"), " Add comparison"),
                   class = "btn-sm btn-outline-primary w-100 mb-1"),
      if (!is.null(reactive_config()))
        actionButton("btn_import_yaml_comps",
                     label = tagList(icon("file-import"), " Import comparisons from YAML"),
                     class = "btn-sm btn-outline-secondary w-100")
      else NULL
    )

    card(
      card_header(icon("wrench"), " Comparisons (manual form)"),
      btns,
      if (length(ids) > 0) tagList(br(), comp_cards) else NULL
    )
  })

  # YAML export button
  output$yaml_export_btn <- renderUI({
    cfg <- reactive_active_config()
    if (is.null(cfg)) return(NULL)
    downloadButton("dl_yaml_export",
                   label = tagList(icon("file-code"), " Export YAML"),
                   class = "btn-outline-info w-100 mt-1 btn-sm")
  })

  # ── OBSERVER: add comparison ─────────────────────────────────────────────────
  observeEvent(input$btn_add_comp, {
    rv$comp_counter <- rv$comp_counter + 1L
    new_id <- rv$comp_counter
    rv$comp_ids <- c(rv$comp_ids, new_id)
    rv$comp_presets[[as.character(new_id)]] <- list()
    log_msg(paste0("Added comparison #", new_id, " to form"), "info")
  })

  # ── OBSERVER: remove comparison ──────────────────────────────────────────────
  observe({
    lapply(rv$comp_ids, function(id) {
      local({
        lid <- id
        observeEvent(input[[paste0("btn_rm_", lid)]], {
          rv$comp_ids    <- setdiff(rv$comp_ids, lid)
          rv$comp_presets[[as.character(lid)]] <- NULL
          log_msg(paste0("Removed comparison #", lid), "info")
        }, ignoreInit = TRUE, ignoreNULL = TRUE, once = TRUE)
      })
    })
  })
  
  # ── OBSERVER: import comparisons from YAML into form ────────────────────────
  observeEvent(input$btn_import_yaml_comps, {
    cfg <- reactive_config()
    req(cfg)
    
    rv$comp_ids     <- integer(0)
    rv$comp_presets <- list()
    
    for (comp in cfg$comparisons) {
      rv$comp_counter <- rv$comp_counter + 1L
      id  <- rv$comp_counter
      rv$comp_ids <- c(rv$comp_ids, id)
      rv$comp_presets[[as.character(id)]] <- comp
    }
    
    updateRadioButtons(session, "config_mode", selected = "form")
    log_msg(paste0("Imported ", length(cfg$comparisons),
                   " comparison(s) from YAML into form"), "ok")
    showNotification(
      paste0(length(cfg$comparisons), " comparison(s) imported from YAML."),
      type = "message", duration = 4
    )
  })
  
  # ── DOWNLOAD: esporta YAML dalla config attiva ───────────────────────────────
  output$dl_yaml_export <- downloadHandler(
    filename = function() {
      cfg <- reactive_active_config()
      paste0(cfg$project_name %||% "config", "_", Sys.Date(), ".yaml")
    },
    content = function(file) {
      cfg <- reactive_active_config()
      req(cfg)
      yaml::write_yaml(cfg, file)
      log_msg(paste0("YAML exported: ", cfg$project_name %||% "—"), "ok")
    }
  )
  
  # ── AUTO-POPOLA FORM quando arriva il YAML ──────────────────────────────────
  observeEvent(reactive_config(), {
    cfg <- reactive_config()
    req(cfg)
    
    updateTextInput(session, "project_name",
                    value = cfg$project_name %||% "")
    
    author_raw <- as.character(cfg$author %||% "")
    if (nchar(author_raw) > 0)
      updateTextInput(session, "author_name", value = author_raw)
    
    gs <- cfg$global_settings
    if (!is.null(gs)) {
      updateNumericInput(session, "undetermined_value",
                         value = gs$undetermined_value %||% 40)
      updateTextInput(session, "housekeeping_genes",
                      value = paste(gs$housekeeping_genes %||% c("GAPDH","HPRT1"),
                                    collapse = ", "))
      updateTextInput(session, "exclude_genes",
                      value = paste(gs$exclude_genes %||% c("RT","gDNA"),
                                    collapse = ", "))
    }
    
    log_msg("Form updated from YAML values", "info")
  })
  
  # ── STATUS BOX ──────────────────────────────────────────────────────────────
  output$status_box <- renderUI({
    data_ok <- !is.null(reactive_data())
    yaml_ok <- !is.null(reactive_config())
    
    data_info <- if (data_ok) {
      df <- reactive_data()
      n_s <- length(unique(df$Sample.Name))
      n_g <- length(unique(df$Target.Name))
      n_u <- sum(df$CT == "Undetermined", na.rm = TRUE)
      paste0(n_s, " samples · ", n_g, " genes · ", n_u, " Undetermined")
    } else ""
    
    yaml_info <- if (yaml_ok) {
      cfg <- reactive_config()
      paste0(length(cfg$comparisons), " comparison(s) defined")
    } else ""
    
    ok_ico   <- tags$span(icon("circle-check"), style = "color:#10b981;")
    off_ico  <- tags$span(icon("circle"), style = "color:#cbd5e1;")
    
    tagList(
      tags$ul(style = "list-style:none; padding-left:0;",
              tags$li(
                if (data_ok)
                  tags$span(ok_ico, " CSV loaded", style = "color:#10b981; font-weight:bold;")
                else
                  tags$span(off_ico, " CSV not loaded", style = "color:#aaa;")
              ),
              if (data_ok) tags$li(tags$small(data_info, style = "color:#555; margin-left:1.6em;")),
              tags$li(
                if (yaml_ok)
                  tags$span(ok_ico, " YAML loaded", style = "color:#10b981; font-weight:bold;")
                else
                  tags$span(off_ico, " YAML not loaded", style = "color:#aaa;")
              ),
              if (yaml_ok) tags$li(tags$small(yaml_info, style = "color:#555; margin-left:1.6em;")),
              tags$li(
                if (length(rv$results) > 0)
                  tags$span(ok_ico, paste0(" ", length(rv$results), " comparison(s) analysed"),
                            style = "color:#10b981; font-weight:bold;")
                else
                  tags$span(off_ico, " No analysis run yet", style = "color:#aaa;")
              )
      )
    )
  })
  
  # ── TAB DATI: summary + tabella ─────────────────────────────────────────────
  output$data_summary <- renderUI({
    df <- reactive_data()
    if (is.null(df)) {
      return(div(class = "alert alert-warning",
                 icon("triangle-exclamation"),
                 " No valid CSV loaded. Please ensure the file has columns: Sample.Name, Target.Name, CT."))
    }
    
    samples  <- sort(unique(df$Sample.Name))
    genes    <- sort(unique(df$Target.Name))
    n_undet  <- sum(df$CT == "Undetermined", na.rm = TRUE)
    ct_num   <- suppressWarnings(as.numeric(gsub(",", ".", df$CT)))
    ct_range <- range(ct_num, na.rm = TRUE)
    
    tagList(
      div(
        class = "gca-stat-row",
        div(class = "gca-stat gca-stat-p",
            div(class = "gca-stat-ico", icon("vials")),
            div(div(class = "gca-stat-val", length(samples)),
                div(class = "gca-stat-lab", "Samples"))),
        div(class = "gca-stat gca-stat-s",
            div(class = "gca-stat-ico", icon("dna")),
            div(div(class = "gca-stat-val", length(genes)),
                div(class = "gca-stat-lab", "Genes"))),
        div(class = "gca-stat gca-stat-w",
            div(class = "gca-stat-ico", icon("circle-xmark")),
            div(div(class = "gca-stat-val", n_undet),
                div(class = "gca-stat-lab", "Undetermined"))),
        div(class = "gca-stat gca-stat-g",
            div(class = "gca-stat-ico", icon("chart-line")),
            div(div(class = "gca-stat-val gca-stat-val-sm",
                    paste0(round(ct_range[1], 1), "–", round(ct_range[2], 1))),
                div(class = "gca-stat-lab", "Range CT")))
      ),
      br(),
      tags$b("Detected samples: "),
      tags$code(paste(samples, collapse = "  |  ")),
      br(), br()
    )
  })
  
  output$table_preview <- renderDT({
    df <- reactive_data()
    req(df)
    cols_show <- intersect(c("Sample.Name","Target.Name","CT","Ct.Mean","Ct.SD"), colnames(df))
    datatable(
      df[, cols_show, drop = FALSE],
      filter  = "top",
      options = list(scrollX = TRUE, pageLength = 15),
      caption = paste0("CT data — ", nrow(df), " total measurements")
    )
  })
  
  # ── TAB CONFIGURAZIONE: summary + cards confronti ───────────────────────────
  output$config_summary <- renderUI({
    cfg <- reactive_config()
    if (is.null(cfg)) {
      return(div(class = "alert alert-info",
                 icon("circle-info"),
                 " No YAML loaded — use the sidebar form to configure comparisons manually."))
    }
    gs <- cfg$global_settings %||% list()
    div(
      div(
        class = "gca-stat-row gca-stat-row-3",
        div(class = "gca-stat gca-stat-p",
            div(class = "gca-stat-ico", icon("flask")),
            div(div(class = "gca-stat-val gca-stat-val-sm", cfg$project_name %||% "—"),
                div(class = "gca-stat-lab", "Project"))),
        div(class = "gca-stat gca-stat-s",
            div(class = "gca-stat-ico", icon("user")),
            div(div(class = "gca-stat-val gca-stat-val-sm", reactive_author()),
                div(class = "gca-stat-lab", "Author"))),
        div(class = "gca-stat gca-stat-g",
            div(class = "gca-stat-ico", icon("layer-group")),
            div(div(class = "gca-stat-val", length(cfg$comparisons)),
                div(class = "gca-stat-lab", "Comparisons")))
      ),
      br(),
      tags$b("Housekeeping: "),
      tags$code(paste(gs$housekeeping_genes %||% "—", collapse = ", ")),
      tags$span("  |  "),
      tags$b("Undetermined value: "),
      tags$code(gs$undetermined_value %||% 40),
      tags$span("  |  "),
      tags$b("Excluded: "),
      tags$code(paste(gs$exclude_genes %||% "—", collapse = ", ")),
      br(), br()
    )
  })
  
  output$comparisons_ui <- renderUI({
    cfg <- reactive_config()
    if (is.null(cfg) || length(cfg$comparisons) == 0) return(NULL)

    # Badge colour metadata for signature types
    sig_meta <- list(
      "all"                  = list(label = "All genes",   color = "#6c757d"),
      "signature_emt"        = list(label = "EMT",          color = "#0d6efd"),
      "signature_metastasis" = list(label = "Metastasis",   color = "#fd7e14"),
      "signature_metabolism" = list(label = "Metabolism",   color = "#198754")
    )

    rows <- lapply(seq_along(cfg$comparisons), function(i) {
      comp     <- cfg$comparisons[[i]]
      sig_type <- comp$signature_type %||% "all"

      # Signature badge with colour
      if (startsWith(sig_type, "custom_")) {
        sig_color <- "#6f42c1"
        sig_lbl   <- paste0("Custom: ", sub("^custom_", "", sig_type))
      } else {
        m         <- sig_meta[[sig_type]] %||% list(label = sig_type, color = "#6c757d")
        sig_color <- m$color; sig_lbl <- m$label
      }
      badge <- tags$span(
        class = "badge",
        style = paste0("background:", sig_color,
                       "; font-size:.70em; padding:3px 8px; border-radius:5px;"),
        sig_lbl
      )

      # Cella campioni
      samp <- comp$samples
      if (is.null(samp) || (length(samp) == 1 && identical(samp, "all"))) {
        samp_cell <- tags$span(tags$em("all"), style = "color:#94a3b8; font-size:.83em;")
      } else {
        samp_cell <- tagList(
          tags$span(paste(samp, collapse = " · "), style = "font-size:.82em;"),
          tags$span(paste0(" (", length(samp), " samp.)"),
                    style = "color:#94a3b8; font-size:.76em;")
        )
      }

      # Cella clustering
      clust_cell <- if (isTRUE(comp$cluster_heatmap))
        tags$span(icon("check"), style = "color:#10b981;", title = "Clustering active")
      else
        tags$span("—", style = "color:#d1d5db;", title = "No clustering")

      tags$tr(
        tags$td(
          tags$span(class = "badge rounded-pill",
                    style = "background:#f1f5f9; color:#64748b; border:1px solid #e2e8f0; font-weight:600;", i),
          style = "width:40px; text-align:center; vertical-align:middle;"
        ),
        tags$td(tags$strong(comp$name, style = "font-size:.88em;")),
        tags$td(badge, style = "vertical-align:middle;"),
        tags$td(
          tags$code(
            style = "font-size:.78em; background:#f1f5f9; color:#0f172a; padding:2px 6px; border-radius:4px; border:1px solid #e2e8f0;",
            comp$reference %||% "—"
          )
        ),
        tags$td(samp_cell, style = "vertical-align:middle;"),
        tags$td(clust_cell, style = "text-align:center; vertical-align:middle; font-size:1em;")
      )
    })

    # Header stile colonne
    th <- function(label, ...) {
      tags$th(label, style = paste0(
        "color:#64748b; font-weight:700; font-size:.72em; ",
        "text-transform:uppercase; letter-spacing:.05em; ",
        "border-bottom:2px solid #e2e8f0; padding:8px 10px;",
        ...
      ))
    }

    card(
      style = "border:1px solid #e3e8f0; border-radius:10px; overflow:hidden;",
      tags$div(
        class = "table-responsive",
        tags$table(
          class = "table table-hover table-sm align-middle mb-0",
          style = "font-size:.88em;",
          tags$thead(
            style = "background:#f8fafc;",
            tags$tr(
              th("#",              "width:40px; text-align:center;"),
              th("Name"),
              th("Gene signature"),
              th("Reference (UT)"),
              th("Samples"),
              th("Cluster",        "width:70px; text-align:center;")
            )
          ),
          tags$tbody(rows)
        )
      )
    )
  })
  
  # ── LOG ─────────────────────────────────────────────────────────────────────
  output$log_output <- renderText({
    if (length(rv$log_lines) == 0) return("No events recorded.")
    paste(rv$log_lines, collapse = "\n")
  })
  
  # ── ANALISI PRINCIPALE ──────────────────────────────────────────────────────
  reactive_analysis <- eventReactive(input$btn_run, {
    df  <- reactive_data()
    cfg <- reactive_active_config()
    
    if (is.null(df))  { showNotification("Please load a valid CSV first.",        type = "warning"); return(list()) }
    if (is.null(cfg)) { showNotification("Configuration not available.",          type = "warning"); return(list()) }
    if (length(cfg$comparisons) == 0) {
      showNotification("No comparisons defined in configuration.", type = "warning")
      return(list())
    }
    
    log_msg(paste0("Analysis started — project: ", cfg$project_name), "info")
    results <- list()
    
    withProgress(message = "Analysis in progress…", value = 0, {
      
      n_comp  <- length(cfg$comparisons)
      ut_val  <- cfg$global_settings$undetermined_value %||% 40
      hk_genes <- cfg$global_settings$housekeeping_genes
      ex_genes <- unique(c(
        cfg$global_settings$exclude_genes %||% character(0),
        c("CTRN", "RQ1", "RQ2", "PCR", "GDNA")
      ))
      
      ct_raw_num <- suppressWarnings(
        as.numeric(gsub(",", ".", df$CT[df$CT != "Undetermined"]))
      )
      ct_max_obs <- suppressWarnings(max(ct_raw_num, na.rm = TRUE))
      if (is.finite(ct_max_obs) && ut_val < ct_max_obs) {
        log_msg(paste0(
          "Warning: undetermined_value (", ut_val,
          ") is lower than the maximum observed CT (", round(ct_max_obs, 1),
          "). High CT values may be treated as expressed."
        ), "warn")
      }
      
      incProgress(0.15, detail = "Normalising CT values…")
      CT_matrix_global <- tryCatch({
        as.data.frame(extract.normalised.CTs(
          df,
          undetermined      = ut_val,
          housekeeping.gene = hk_genes,
          exclude_genes     = ex_genes
        ))
      }, error = function(e) {
        log_msg(paste0("CT normalisation error: ", e$message), "error")
        NULL
      })
      
      if (is.null(CT_matrix_global)) {
        showNotification("CT normalisation failed. Check the Log tab.",
                         type = "error", duration = 8)
        return(list())
      }
      
      log_msg(paste0("Normalised CT matrix: ",
                     nrow(CT_matrix_global), " genes × ",
                     ncol(CT_matrix_global), " samples"), "ok")
      
      all_signatures <- get_signatures(cfg)
      
      for (i in seq_along(cfg$comparisons)) {
        comp <- cfg$comparisons[[i]]
        incProgress(0.85 / n_comp, detail = paste("Processing:", comp$name))
        
        CT_matrix <- CT_matrix_global
        
        if (length(comp$samples) == 1 && comp$samples == "all") {
          selected_matrix <- CT_matrix
        } else {
          avail <- intersect(comp$samples, colnames(CT_matrix))
          if (length(avail) == 0) {
            log_msg(paste0("[", comp$name, "] No samples found: ",
                           paste(comp$samples, collapse = ", ")), "warn")
            next
          }
          if (length(avail) < length(comp$samples)) {
            missing_s <- setdiff(comp$samples, colnames(CT_matrix))
            log_msg(paste0("[", comp$name, "] Missing samples: ",
                           paste(missing_s, collapse = ", ")), "warn")
          }
          selected_matrix <- CT_matrix[, avail, drop = FALSE]
        }
        
        sig_genes <- if (comp$signature_type == "all") {
          rownames(selected_matrix)
        } else if (!is.null(all_signatures[[comp$signature_type]])) {
          all_signatures[[comp$signature_type]]
        } else {
          log_msg(paste0("[", comp$name, "] Signature '", comp$signature_type,
                         "' not found — using all genes"), "warn")
          rownames(selected_matrix)
        }
        
        genes_to_keep <- setdiff(
          intersect(sig_genes, rownames(selected_matrix)),
          cfg$global_settings$exclude_genes %||% character(0)
        )
        
        if (length(genes_to_keep) == 0) {
          log_msg(paste0("[", comp$name, "] No genes remaining after filters — skipping"), "warn")
          next
        }
        
        final_matrix <- selected_matrix[genes_to_keep, , drop = FALSE]
        final_matrix <- final_matrix[order(rownames(final_matrix)), , drop = FALSE]
        
        ut_ref    <- comp$reference
        log2fc_df <- NULL
        
        if (ut_ref %in% colnames(final_matrix)) {
          log2fc_mat <- log2(as.matrix(final_matrix) / final_matrix[, ut_ref])
          log2fc_mat[is.infinite(log2fc_mat)] <- NA
          log2fc_df <- as.data.frame(log2fc_mat) %>%
            rownames_to_column("Gene") %>%
            arrange(Gene)
        } else {
          log_msg(paste0("[", comp$name, "] Reference '", ut_ref,
                         "' not found — log2FC not computed"), "warn")
        }
        
        # Z-score su scala log2 (= −deltaCT): log2(2^−deltaCT) = −deltaCT
        # Avoid z-scoring raw 2^(-deltaCT) values which are on an exponential scale
        mat_log2      <- log2(as.matrix(final_matrix))
        mat_log2[!is.finite(mat_log2)] <- NA
        matrix_scaled <- t(scale(t(mat_log2)))
        
        results[[comp$name]] <- list(
          matrix        = final_matrix,
          matrix_scaled = matrix_scaled,
          log2fc_df     = log2fc_df,
          info          = comp
        )
        
        log_msg(paste0("[", comp$name, "] OK — ",
                       nrow(final_matrix), " genes · ",
                       ncol(final_matrix), " samples"), "ok")
      }
    })
    
    if (length(results) == 0)
      showNotification("No comparisons completed. Check the Log tab.", type = "warning")
    else
      showNotification(paste0("Analysis complete: ", length(results), " comparison(s)"),
                       type = "message", duration = 5)
    
    results
  })
  
  observe({
    res <- reactive_analysis()
    rv$results <- res
  })
  
  observeEvent(input$btn_run, {
    df  <- reactive_data()
    cfg <- reactive_active_config()
    if (is.null(df))  { showNotification("Please load a valid CSV first.",   type = "warning"); return() }
    if (is.null(cfg)) { showNotification("Configuration not available.",     type = "warning"); return() }
  }, ignoreInit = TRUE)
  
  # ===========================================================================
  # RISULTATI — Opzione B: dropdown confronto + UNA navbar contenuti
  # ===========================================================================
  output$results_ui <- renderUI({
    res <- rv$results
    if (length(res) == 0) {
      return(card(
        p(style = "color:grey; padding:1.5em;",
          icon("circle-info"), " No results available. ",
          "Load your data and press ", tags$b("Run Analysis"), ".")
      ))
    }
    
    comp_names <- names(res)
    
    tagList(
      div(
        class = "gca-results-toolbar",
        div(
          class = "gca-comp-picker",
          tags$label("Comparison", class = "gca-picker-label"),
          selectInput("sel_comp", label = NULL,
                      choices  = comp_names,
                      selected = comp_names[1],
                      width    = "260px")
        ),
        div(class = "gca-toolbar-sep"),
        div(
          class = "gca-results-dl",
          downloadButton("dl_png_sel",      tagList(icon("image"),      " PNG"),
                         class = "btn-sm btn-outline-primary"),
          downloadButton("dl_pdf_sel",      tagList(icon("file-pdf"),   " PDF z-score"),
                         class = "btn-sm btn-outline-secondary"),
          downloadButton("dl_pdf_norm_sel", tagList(icon("file-pdf"),   " PDF norm."),
                         class = "btn-sm btn-outline-secondary"),
          downloadButton("dl_excel_sel",    tagList(icon("file-excel"), " Excel"),
                         class = "btn-sm btn-outline-success"),
          downloadButton("dl_txt_sel",      tagList(icon("file-lines"), " Matrice .txt"),
                         class = "btn-sm btn-outline-dark"),
          downloadButton("dl_pdf_log2fc_sel", tagList(icon("file-pdf"), " PDF log2FC"),
                         class = "btn-sm btn-outline-secondary")
        )
      ),
      
      div(
        class = "gca-results-vizopts",
        # ── Clustering ──────────────────────────────────────────────────────────
        div(
          class = "gca-vopts-group",
          tags$span(class = "gca-picker-label", icon("sitemap"), " Clustering:"),
          radioButtons("ht_cluster", label = NULL,
                       choices  = c("None"        = "none",
                                    "Rows"        = "rows",
                                    "Columns"     = "cols",
                                    "Rows+Cols"   = "both"),
                       selected = "none",
                       inline   = TRUE)
        ),
        div(class = "gca-vopts-sep"),
        # ── Ordine colonne ──────────────────────────────────────────────────────
        div(
          class = "gca-vopts-group",
          tags$span(class = "gca-picker-label", icon("arrow-right-arrow-left"), " Column order:"),
          uiOutput("col_order_sel_ui")
        )
      ),

      uiOutput("results_statrow"),

      navset_card_underline(
        id          = "results_content_tabs",
        full_screen = TRUE,
        nav_panel(
          title = "Tabella log2FC",
          icon  = icon("table"),
          br(),
          uiOutput("results_table_or_warn")
        ),
        nav_panel(
          title = "Heatmap Z-score",
          icon  = icon("map"),
          br(),
          plotOutput("ht_zscore_sel", height = "auto")
        ),
        nav_panel(
          title = "Normalised Heatmap",
          icon  = icon("layer-group"),
          br(),
          plotOutput("ht_norm_sel", height = "auto")
        ),
        nav_panel(
          title = "Heatmap Log2FC",
          icon  = icon("chart-bar"),
          br(),
          uiOutput("results_log2fc_warn_or_plot")
        )
      )
    )
  })
  
  current_result <- reactive({
    res <- rv$results
    req(length(res) > 0)
    sel <- input$sel_comp %||% names(res)[1]
    if (!sel %in% names(res)) sel <- names(res)[1]
    res[[sel]]
  })

  output$col_order_sel_ui <- renderUI({
    r    <- current_result()
    samp <- colnames(r$matrix)
    selectizeInput("ht_col_order", label = NULL,
                   choices  = samp,
                   selected = samp,
                   multiple = TRUE,
                   width    = "360px",
                   options  = list(
                     placeholder = "Select in desired order…",
                     plugins     = list("remove_button")
                   ))
  })

  output$results_statrow <- renderUI({
    r <- current_result()
    div(
      class = "gca-stat-row",
      div(class = "gca-stat gca-stat-p",
          div(class = "gca-stat-ico", icon("dna")),
          div(div(class = "gca-stat-val", nrow(r$matrix)),
              div(class = "gca-stat-lab", "Genes analysed"))),
      div(class = "gca-stat gca-stat-s",
          div(class = "gca-stat-ico", icon("vials")),
          div(div(class = "gca-stat-val", ncol(r$matrix)),
              div(class = "gca-stat-lab", "Samples"))),
      div(class = "gca-stat gca-stat-i",
          div(class = "gca-stat-ico", icon("tag")),
          div(div(class = "gca-stat-val gca-stat-val-sm",
                  r$info$signature_type %||% "all"),
              div(class = "gca-stat-lab", "Signature"))),
      div(class = "gca-stat gca-stat-g",
          div(class = "gca-stat-ico", icon("anchor")),
          div(div(class = "gca-stat-val gca-stat-val-sm",
                  r$info$reference %||% "—"),
              div(class = "gca-stat-lab", "Reference")))
    )
  })
  
  output$results_table_or_warn <- renderUI({
    r <- current_result()
    if (!is.null(r$log2fc_df)) {
      DTOutput("dt_log2fc_sel")
    } else {
      div(class = "alert alert-warning",
          icon("triangle-exclamation"),
          " log2FC not computed: reference sample '",
          r$info$reference %||% "—",
          "' not found in the matrix.")
    }
  })
  
  # ── RENDER + DOWNLOAD for selected comparison (fixed IDs *_sel) ─────────────

  # Colori reattivi dalle preferenze utente (sidebar)
  reactive_ht_colors <- reactive({
    c(
      input$col_low  %||% "#4575b4",
      input$col_mid  %||% "#f7f7f7",
      input$col_high %||% "#d73027"
    )
  })

  .norm_col_fun <- function(mat_nr, cols = c("#4575b4", "#f7f7f7", "#d73027")) {
    norm_vec <- as.vector(mat_nr)
    norm_vec <- norm_vec[is.finite(norm_vec) & norm_vec > 0]
    if (length(norm_vec) >= 3) {
      colorRamp2(quantile(norm_vec, c(0.1, 0.5, 0.9)), cols)
    } else {
      colorRamp2(c(0, 0.5, 1), cols)
    }
  }
  .zs_col_fun <- function(cols = c("#4575b4", "#f7f7f7", "#d73027")) {
    colorRamp2(c(-2, 0, 2), cols)
  }
  .lfc_col_fun <- function(mat_lfc, cols = c("#4575b4", "#f7f7f7", "#d73027")) {
    vals <- mat_lfc[is.finite(mat_lfc)]
    lim  <- if (length(vals) > 0) max(2, ceiling(max(abs(vals), na.rm = TRUE))) else 2
    colorRamp2(c(-lim, 0, lim), cols)
  }

  # Reorder matrix columns according to user input
  .reorder_cols <- function(mat, col_order) {
    if (!is.null(col_order) && length(col_order) > 0) {
      valid <- intersect(col_order, colnames(mat))
      if (length(valid) == ncol(mat)) return(mat[, valid, drop = FALSE])
    }
    mat
  }
  
  output$dt_log2fc_sel <- renderDT({
    r <- current_result()
    req(!is.null(r$log2fc_df))
    df_show  <- r$log2fc_df
    num_cols <- setdiff(colnames(df_show), "Gene")
    datatable(
      df_show,
      rownames = FALSE,
      filter   = "top",
      options  = list(scrollX = TRUE, pageLength = 15,
                      order = list(list(1, "desc")))
    ) %>%
      formatRound(columns = num_cols, digits = 3) %>%
      formatStyle(
        columns = num_cols,
        backgroundColor = styleInterval(
          cuts   = c(-1, -0.5, 0, 0.5, 1),
          values = c("#4a90d9","#a8c8f0","#f5f5f5","#f5f5f5","#f8b4a0","#e05c3a")
        )
      )
  })
  
  output$ht_zscore_sel <- renderPlot({
    r      <- current_result()
    cfg    <- reactive_active_config()
    author <- reactive_author()
    cols   <- reactive_ht_colors()
    mat    <- .reorder_cols(r$matrix_scaled, input$ht_col_order)
    draw_heatmap_cht(mat, .zs_col_fun(cols), "z-score", "Z-score",
                     input$sel_comp, r$info, cfg, author,
                     cluster_override = input$ht_cluster %||% "none")
  }, height = function() {
    r <- current_result()
    plot_h_px(nrow(r$matrix), r$info$signature_type %||% "all")
  })
  
  output$ht_norm_sel <- renderPlot({
    r      <- current_result()
    cfg    <- reactive_active_config()
    author <- reactive_author()
    cols   <- reactive_ht_colors()
    mat_nr <- .reorder_cols(as.matrix(r$matrix), input$ht_col_order)
    draw_heatmap_cht(mat_nr, .norm_col_fun(mat_nr, cols), "Norm.", "Norm. CT",
                     paste0(input$sel_comp, " — Normalized"), r$info, cfg, author,
                     cluster_override = input$ht_cluster %||% "none")
  }, height = function() {
    r <- current_result()
    plot_h_px(nrow(r$matrix), r$info$signature_type %||% "all")
  })
  
  output$results_log2fc_warn_or_plot <- renderUI({
    r <- current_result()
    if (!is.null(r$log2fc_df)) {
      plotOutput("ht_log2fc_sel", height = "auto")
    } else {
      div(class = "alert alert-warning m-3",
          icon("triangle-exclamation"),
          " Log2FC not computed: reference sample '",
          r$info$reference %||% "—",
          "' was not found in the matrix.")
    }
  })

  output$ht_log2fc_sel <- renderPlot({
    r      <- current_result()
    req(!is.null(r$log2fc_df))
    cfg    <- reactive_active_config()
    author <- reactive_author()
    cols      <- reactive_ht_colors()
    log2fc_df <- r$log2fc_df
    mat_lfc   <- as.matrix(log2fc_df[, -1, drop = FALSE])
    rownames(mat_lfc) <- log2fc_df$Gene
    if (ncol(mat_lfc) <= 2) {
      # 2 condizioni: bar chart ranked — molto più leggibile della heatmap
      print(draw_log2fc_bars(log2fc_df, r$info, cfg, author))
    } else {
      mat_lfc <- .reorder_cols(mat_lfc, input$ht_col_order)
      draw_heatmap_cht(mat_lfc, .lfc_col_fun(mat_lfc, cols), "log2FC", "Log2 Fold Change",
                       paste0(input$sel_comp, " — Log2FC"), r$info, cfg, author,
                       cluster_override = input$ht_cluster %||% "none")
    }
  }, height = function() {
    r <- current_result()
    plot_h_px(nrow(r$matrix), r$info$signature_type %||% "all")
  })

  # Prevent Shiny from suspending these outputs when their tab is inactive.
  # Without this, plots in non-active nav_panels never render on first switch.
  outputOptions(output, "ht_zscore_sel",  suspendWhenHidden = FALSE)
  outputOptions(output, "ht_norm_sel",    suspendWhenHidden = FALSE)
  outputOptions(output, "ht_log2fc_sel",  suspendWhenHidden = FALSE)
  outputOptions(output, "dt_log2fc_sel",  suspendWhenHidden = FALSE)

  output$dl_png_sel <- downloadHandler(
    filename = function() paste0("Heatmap_", input$sel_comp, "_zscore_", Sys.Date(), ".png"),
    content  = function(file) {
      r     <- current_result(); cfg <- reactive_active_config(); author <- reactive_author()
      cols  <- reactive_ht_colors()
      sig_t <- r$info$signature_type %||% "all"
      mat   <- .reorder_cols(r$matrix_scaled, input$ht_col_order)
      h_in  <- plot_h_px(nrow(mat), sig_t, dpi = 300) / 300
      png(file, width = 14 * 300, height = h_in * 300, res = 300, type = "cairo")
      showtext_begin()
      draw_heatmap_cht(mat, .zs_col_fun(cols), "z-score", "Z-score",
                       input$sel_comp, r$info, cfg, author, font_mult = 3,
                       cluster_override = input$ht_cluster %||% "none")
      showtext_end(); dev.off()
    }
  )
  
  output$dl_pdf_sel <- downloadHandler(
    filename = function() paste0("Heatmap_", input$sel_comp, "_zscore_", Sys.Date(), ".pdf"),
    content  = function(file) {
      r    <- current_result(); cfg <- reactive_active_config(); author <- reactive_author()
      cols <- reactive_ht_colors()
      sig_t <- r$info$signature_type %||% "all"
      mat  <- .reorder_cols(r$matrix_scaled, input$ht_col_order)
      h_in <- plot_h_px(nrow(mat), sig_t, dpi = 96) / 96
      cairo_pdf(file, width = 14, height = h_in)
      showtext_begin()
      draw_heatmap_cht(mat, .zs_col_fun(cols), "z-score", "Z-score",
                       input$sel_comp, r$info, cfg, author,
                       cluster_override = input$ht_cluster %||% "none")
      showtext_end(); dev.off()
    }
  )
  
  output$dl_pdf_norm_sel <- downloadHandler(
    filename = function() paste0("Heatmap_", input$sel_comp, "_norm_", Sys.Date(), ".pdf"),
    content  = function(file) {
      r    <- current_result(); cfg <- reactive_active_config(); author <- reactive_author()
      cols <- reactive_ht_colors()
      sig_t  <- r$info$signature_type %||% "all"
      mat_nr <- .reorder_cols(as.matrix(r$matrix), input$ht_col_order)
      h_in   <- plot_h_px(nrow(mat_nr), sig_t, dpi = 96) / 96
      cairo_pdf(file, width = 14, height = h_in)
      showtext_begin()
      draw_heatmap_cht(mat_nr, .norm_col_fun(mat_nr, cols), "Norm.", "Norm. CT",
                       paste0(input$sel_comp, " — Normalized"), r$info, cfg, author,
                       cluster_override = input$ht_cluster %||% "none")
      showtext_end(); dev.off()
    }
  )
  
  output$dl_excel_sel <- downloadHandler(
    filename = function() paste0("log2FC_", input$sel_comp, "_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      r <- current_result(); req(!is.null(r$log2fc_df))
      cfg <- reactive_active_config(); gs <- cfg$global_settings %||% list()
      legend_df <- data.frame(
        Parameter     = c("log2FC > 0", "log2FC < 0", "log2FC = 0",
                          "Interpretation", "Reference", "Housekeeping"),
        Value         = c("Upregulated", "Downregulated", "No Change",
                          "Each unit = 2-fold difference (e.g. log2FC 2 = 4x more expressed)",
                          r$info$reference %||% "—",
                          paste(gs$housekeeping_genes %||% "—", collapse = ", "))
      )
      write_xlsx(list(Log2FC_Data = r$log2fc_df, Legend = legend_df), file)
    }
  )
  
  output$dl_txt_sel <- downloadHandler(
    filename = function() paste0("Scaled_Matrix_", input$sel_comp, "_", Sys.Date(), ".txt"),
    content  = function(file) {
      r <- current_result()
      write.table(r$matrix_scaled, file, sep = "\t", quote = FALSE)
    }
  )

  output$dl_pdf_log2fc_sel <- downloadHandler(
    filename = function() paste0("Log2FC_", input$sel_comp, "_", Sys.Date(), ".pdf"),
    content  = function(file) {
      r <- current_result(); req(!is.null(r$log2fc_df))
      cfg       <- reactive_active_config(); author <- reactive_author()
      cols      <- reactive_ht_colors()
      sig_t     <- r$info$signature_type %||% "all"
      log2fc_df <- r$log2fc_df
      mat_lfc   <- as.matrix(log2fc_df[, -1, drop = FALSE])
      rownames(mat_lfc) <- log2fc_df$Gene
      h_in <- plot_h_px(nrow(mat_lfc), sig_t, dpi = 96) / 96
      cairo_pdf(file, width = 10, height = h_in)
      showtext_begin()
      if (ncol(mat_lfc) <= 2) {
        print(draw_log2fc_bars(log2fc_df, r$info, cfg, author))
      } else {
        mat_lfc <- .reorder_cols(mat_lfc, input$ht_col_order)
        draw_heatmap_cht(mat_lfc, .lfc_col_fun(mat_lfc, cols), "log2FC", "Log2 Fold Change",
                         paste0(input$sel_comp, " — Log2FC"), r$info, cfg, author,
                         cluster_override = input$ht_cluster %||% "none")
      }
      showtext_end(); dev.off()
    }
  )

  # ── TAB DOWNLOAD: ZIP globale + riepilogo ───────────────────────────────────
  output$downloads_ui <- renderUI({
    res <- rv$results
    if (length(res) == 0) {
      return(card(
        p(style = "color:grey; padding:1.5em;",
          icon("circle-info"),
          " No files available. Please run the analysis first.")
      ))
    }
    
    header_card <- card(
      card_header(
        class = "d-flex justify-content-between align-items-center",
        tags$span(icon("box-archive"), " ", tags$b("Full Download")),
        downloadButton("dl_zip_all",
                       label = tagList(icon("file-zipper"), " Download all (ZIP)"),
                       class = "btn-success")
      ),
      p(style = "color:#555;",
        "The ZIP contains, for each comparison: heatmap PNG (300 DPI), PDF z-score, ",
        "normalised PDF, Excel log2FC table, and scaled matrix (.txt). ",
        "To download individual files for a comparison, use the buttons in the ",
        tags$b("Results"), " tab after selecting the comparison.")
    )
    
    summary_rows <- lapply(names(res), function(cn) {
      r <- res[[cn]]
      tags$tr(
        tags$td(tags$b(cn)),
        tags$td(paste0(nrow(r$matrix), " genes")),
        tags$td(paste0(ncol(r$matrix), " samples")),
        tags$td(r$info$signature_type %||% "all"),
        tags$td(
          if (!is.null(r$log2fc_df))
            tags$span(icon("circle-check"), " yes", style = "color:#10b981;")
          else
            tags$span(icon("circle-xmark"), " no", style = "color:#aaa;")
        )
      )
    })
    
    detail_card <- card(
      card_header(icon("list"), " Comparison summary"),
      tags$table(
        class = "table table-sm table-hover mb-0",
        tags$thead(
          tags$tr(
            tags$th("Comparison"), tags$th("Genes"), tags$th("Samples"),
            tags$th("Signature"), tags$th("log2FC")
          )
        ),
        tags$tbody(summary_rows)
      )
    )
    
    tagList(header_card, br(), detail_card)
  })
  
  # ── DOWNLOAD ZIP (tutti i file) ──────────────────────────────────────────────
  output$dl_zip_all <- downloadHandler(
    filename = function() {
      cfg <- reactive_active_config()
      paste0(cfg$project_name %||% "GeneCard", "_", Sys.Date(), ".zip")
    },
    content = function(zip_file) {
      res    <- rv$results
      cfg    <- reactive_active_config()
      author <- reactive_author()
      
      tmpdir <- file.path(tempdir(),
                          paste0("gca_", format(Sys.time(), "%H%M%S")))
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
      
      withProgress(message = "Generating ZIP…", value = 0, {
        n <- length(res)
        
        for (i in seq_along(res)) {
          cn <- names(res)[i]
          r  <- res[[cn]]
          incProgress(1 / n, detail = cn)
          
          comp_dir <- file.path(tmpdir, cn)
          dir.create(comp_dir, showWarnings = FALSE)
          
          if (!is.null(r$log2fc_df)) {
            gs <- cfg$global_settings %||% list()
            legend_df <- data.frame(
              Parameter     = c("log2FC > 0", "log2FC < 0", "log2FC = 0",
                                "Interpretation", "Reference", "Housekeeping"),
              Value         = c("Upregulated", "Downregulated", "No Change",
                                "Each unit = 2-fold difference",
                                r$info$reference %||% "—",
                                paste(gs$housekeeping_genes %||% "—", collapse = ", "))
            )
            write_xlsx(
              list(Log2FC_Data = r$log2fc_df, Legend = legend_df),
              file.path(comp_dir, paste0("log2FC_", cn, ".xlsx"))
            )
          }
          
          write.table(r$matrix_scaled,
                      file.path(comp_dir, paste0("Scaled_Matrix_", cn, ".txt")),
                      sep = "\t", quote = FALSE)
          
          col_zs <- colorRamp2(c(-2, 0, 2), c("dodgerblue4", "gray95", "firebrick3"))
          h_in   <- plot_h_px(nrow(r$matrix_scaled),
                              r$info$signature_type %||% "all", dpi = 300) / 300
          png_path <- file.path(comp_dir, paste0("Heatmap_", cn, "_zscore.png"))
          png(png_path, width = 14 * 300, height = h_in * 300, res = 300, type = "cairo")
          showtext_begin()
          draw_heatmap_cht(r$matrix_scaled, col_zs, "z-score", "Z-score",
                           cn, r$info, cfg, author, font_mult = 3)
          showtext_end()
          dev.off()
          
          h_in_pdf <- plot_h_px(nrow(r$matrix_scaled),
                                r$info$signature_type %||% "all", dpi = 96) / 96
          pdf_path <- file.path(comp_dir, paste0("Heatmap_", cn, "_zscore.pdf"))
          cairo_pdf(pdf_path, width = 14, height = h_in_pdf)
          showtext_begin()
          draw_heatmap_cht(r$matrix_scaled, col_zs, "z-score", "Z-score",
                           cn, r$info, cfg, author)
          showtext_end()
          dev.off()
          
          mat_nr   <- as.matrix(r$matrix)
          norm_vec <- as.vector(mat_nr)
          norm_vec <- norm_vec[is.finite(norm_vec) & norm_vec > 0]
          col_nr <- if (length(norm_vec) >= 3) {
            colorRamp2(
              quantile(norm_vec, c(0.1, 0.5, 0.9)),
              c("dodgerblue4", "gray95", "firebrick3")
            )
          } else {
            colorRamp2(c(0, 0.5, 1), c("dodgerblue4", "gray95", "firebrick3"))
          }
          pdf_norm <- file.path(comp_dir, paste0("Heatmap_", cn, "_norm.pdf"))
          cairo_pdf(pdf_norm, width = 14, height = h_in_pdf)
          showtext_begin()
          draw_heatmap_cht(mat_nr, col_nr, "Norm.", "Norm. CT",
                           paste0(cn, " — Normalized"), r$info, cfg, author)
          showtext_end()
          dev.off()
          
          log_msg(paste0("[ZIP] ", cn, " — files generated"), "ok")
        }
      })
      
      owd <- setwd(tmpdir)
      on.exit(setwd(owd), add = TRUE)
      all_files <- list.files(".", recursive = TRUE, full.names = FALSE)
      utils::zip(zipfile = zip_file, files = all_files)
      
      log_msg(paste0("ZIP created: ", length(all_files), " files"), "ok")
    }
  )
  
}

# ==============================================================================
shinyApp(ui = ui, server = server)