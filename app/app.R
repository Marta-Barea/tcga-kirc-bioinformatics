suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(survival)
})

format_metric <- function(value) {
  format(value, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
}

get_variable_label <- function(variable_name, labels) {
  if (variable_name %in% names(labels)) {
    return(labels[[variable_name]])
  }

  gsub("_", " ", variable_name)
}

resolve_project_dir <- function() {
  app_dir <- shiny::getShinyOption("appDir")

  if (!is.null(app_dir) && nzchar(app_dir)) {
    return(dirname(normalizePath(app_dir, winslash = "/", mustWork = FALSE)))
  }

  frame_files <- vapply(
    sys.frames(),
    function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else as.character(ofile)
    },
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]

  if (length(frame_files) > 0) {
    return(dirname(dirname(normalizePath(frame_files[[length(frame_files)]], winslash = "/", mustWork = FALSE))))
  }

  project_dir <- getwd()

  if (basename(project_dir) %in% c("scripts", "app")) {
    project_dir <- dirname(project_dir)
  }

  project_dir
}

build_clinical_data <- function(clinical_path) {
  clinical_raw <- readRDS(clinical_path)
  age_median <- median(clinical_raw$age_at_diagnosis, na.rm = TRUE)
  days_per_year <- 365.25
  days_per_month <- days_per_year / 12

  clinical_raw %>%
    transmute(
      patient_id = submitter_id,
      age_at_diagnosis = dplyr::coalesce(age_at_diagnosis, age_median) / days_per_year,
      ajcc_pathologic_stage = if_else(
        is.na(ajcc_pathologic_stage),
        "Unknown",
        ajcc_pathologic_stage
      ),
      tumor_grade = if_else(
        is.na(tumor_grade),
        "GX",
        tumor_grade
      ),
      sex = dplyr::coalesce(sex_at_birth, gender, "Unknown"),
      vital_status = dplyr::coalesce(vital_status, "Unknown"),
      death_event = if_else(is.na(days_to_death), 0L, 1L),
      survival_time_months = if_else(
        is.na(days_to_death),
        days_to_last_follow_up,
        days_to_death
      ) / days_per_month
    ) %>%
    mutate(
      survival_time_months = if_else(
        is.na(survival_time_months),
        0,
        survival_time_months
      )
    )
}

build_expression_data <- function(se_path) {
  se_object <- readRDS(se_path)
  assay_name <- if ("tpm_unstrand" %in% assayNames(se_object)) {
    "tpm_unstrand"
  } else {
    assayNames(se_object)[1]
  }

  expression_matrix <- assay(se_object, assay_name)
  sample_ids <- colnames(expression_matrix)
  gene_names <- rowData(se_object)$gene_name

  if (is.null(gene_names)) {
    gene_names <- rownames(expression_matrix)
  }

  gene_names <- ifelse(is.na(gene_names) | gene_names == "", rownames(expression_matrix), gene_names)
  unique_gene_names <- make.unique(as.character(gene_names))
  rownames(expression_matrix) <- unique_gene_names

  sample_metadata <- tibble::tibble(
    sample_id = sample_ids,
    patient_id = substr(sample_ids, 1, 12),
    sample_type_code = substr(sample_ids, 14, 15),
    tissue_group = dplyr::case_when(
      sample_type_code == "01" ~ "Tumor primario",
      sample_type_code == "11" ~ "Tejido normal",
      TRUE ~ "Otro"
    )
  ) %>%
    filter(sample_type_code %in% c("01", "11"))

  expression_matrix <- expression_matrix[, sample_metadata$sample_id, drop = FALSE]

  list(
    matrix = expression_matrix,
    metadata = sample_metadata,
    assay_name = assay_name,
    gene_choices = sort(unique(rownames(expression_matrix)))
  )
}

build_survival_dataset <- function(clinical_data, sample_data, expression_matrix, gene_symbol) {
  tumor_samples <- sample_data %>%
    filter(tissue_group == "Tumor primario")

  if (!(gene_symbol %in% rownames(expression_matrix))) {
    return(clinical_data %>% filter(patient_id %in% tumor_samples$patient_id))
  }

  gene_values <- tibble::tibble(
    patient_id = tumor_samples$patient_id,
    gene_expression = as.numeric(expression_matrix[gene_symbol, tumor_samples$sample_id]),
    sample_id = tumor_samples$sample_id
  ) %>%
    group_by(patient_id) %>%
    summarise(gene_expression = mean(gene_expression, na.rm = TRUE), .groups = "drop")

  survival_data <- clinical_data %>%
    filter(patient_id %in% tumor_samples$patient_id) %>%
    inner_join(gene_values, by = "patient_id")

  gene_cutoff <- median(survival_data$gene_expression, na.rm = TRUE)

  survival_data %>%
    mutate(
      gene_group = if_else(gene_expression >= gene_cutoff, "Alta expresión", "Baja expresión")
    )
}

make_km_data <- function(fit_object) {
  fit_summary <- summary(fit_object)

  if (length(fit_summary$time) == 0) {
    return(tibble::tibble())
  }

  strata_names <- fit_summary$strata

  if (is.null(strata_names)) {
    strata_names <- rep("Cohorte completa", length(fit_summary$time))
  }

  tibble::tibble(
    time = fit_summary$time,
    surv = fit_summary$surv,
    strata = gsub(".*=", "", strata_names)
  ) %>%
    group_by(strata) %>%
    arrange(time, .by_group = TRUE) %>%
    reframe(
      time = c(0, time),
      surv = c(1, surv)
    )
}

project_dir <- resolve_project_dir()
clinical_path <- file.path(project_dir, "data", "prepared", "clinical_data_TCGA_KIRC.rds")
se_path <- file.path(project_dir, "data", "prepared", "tcga_kirc_star_counts_se.rds")

if (!file.exists(clinical_path) || !file.exists(se_path)) {
  stop(
    "No se encontraron los datos preparados. Ejecuta primero scripts/download_tcga_kirc_data.R.",
    call. = FALSE
  )
}

clinical_data <- build_clinical_data(clinical_path)
expression_bundle <- build_expression_data(se_path)
sample_data <- expression_bundle$metadata
expression_matrix <- expression_bundle$matrix
gene_choices <- expression_bundle$gene_choices
default_gene <- if ("CA9" %in% gene_choices) "CA9" else gene_choices[1]

analysis_data <- sample_data %>%
  left_join(clinical_data, by = "patient_id")

candidate_genes <- c(
  "CA9",
  "NDUFA4L2",
  "EGLN3",
  "HILPDA",
  "SCARB1",
  "STC2",
  "COL23A1",
  "CDKN2A",
  "GABRD",
  "PVT1",
  "TTC21B-AS1",
  "LINC00887"
)

gene_choices <- candidate_genes[candidate_genes %in% gene_choices]
if (length(gene_choices) == 0) {
  stop("No se encontraron en la matriz de expresión los genes candidatos del informe.", call. = FALSE)
}

theme_set(
  theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, colour = "#163a5f"),
      plot.subtitle = element_text(colour = "#5b6f82"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#d9e7f5", linewidth = 0.45),
      panel.grid.major.y = element_line(colour = "#edf4fb", linewidth = 0.45),
      plot.background = element_rect(fill = "#ffffff", colour = NA),
      panel.background = element_rect(fill = "#ffffff", colour = NA),
      legend.position = "bottom"
    )
)

app_theme <- bs_theme(
  version = 5,
  bg = "#ffffff",
  fg = "#16324f",
  primary = "#1f5f99",
  secondary = "#6c84a1",
  success = "#238a8d",
  info = "#3b82c4",
  base_font = font_collection("Avenir Next", "Helvetica Neue", "Segoe UI", "Arial", "sans-serif"),
  heading_font = font_collection("Avenir Next", "Helvetica Neue", "Segoe UI", "Arial", "sans-serif")
)

ui <- page_navbar(
  title = "TCGA-KIRC Explorer",
  theme = app_theme,
  fillable = FALSE,
  header = tagList(
    tags$head(
      tags$style(HTML(
      ":root {
         --app-bg: #ffffff;
         --surface: #f8fbff;
         --surface-strong: #eef6ff;
         --border-soft: rgba(46, 93, 140, 0.12);
         --ink: #16324f;
         --muted: #5b6f82;
         --blue: #1f5f99;
         --blue-deep: #163a5f;
         --viridis: #238a8d;
         --viridis-mid: #2d708e;
         --viridis-deep: #355f8d;
         --viridis-bright: #1f9e89;
       }
       body {
         background: linear-gradient(180deg, #fbfdff 0%, #ffffff 45%);
         color: var(--ink);
       }
       .navbar {
         background: rgba(255,255,255,0.92) !important;
         backdrop-filter: blur(8px);
         border-bottom: 1px solid var(--border-soft);
       }
       .navbar-default .navbar-nav > li > a,
       .navbar-default .navbar-brand {
         color: var(--blue-deep) !important;
       }
       .navbar-default .navbar-nav > .active > a,
       .navbar-default .navbar-nav > .active > a:focus,
       .navbar-default .navbar-nav > .active > a:hover {
         color: var(--blue) !important;
         background: transparent !important;
         box-shadow: inset 0 -2px 0 var(--blue);
       }
       .hero-banner {
         position: relative;
         margin: 1rem 0 1.75rem 0;
         padding: 1.35rem 0 1.1rem 0;
         background: transparent;
         overflow: hidden;
       }
       .hero-banner::before {
         display: none;
         content: none;
       }
       .hero-banner::after {
         display: none;
         content: none;
       }
       .hero-copy span {
         position: relative;
         display: inline-flex;
         align-items: center;
         gap: 0.55rem;
         text-transform: uppercase;
         letter-spacing: 0.12em;
         font-size: 0.78rem;
         color: var(--viridis-mid);
         font-weight: 700;
       }
       .hero-copy span::before {
         content: '';
         width: 2.4rem;
         height: 2px;
         background: linear-gradient(90deg, var(--viridis-deep), var(--viridis-bright));
         border-radius: 999px;
       }
       .hero-copy h1 {
         position: relative;
         margin-top: 0.25rem;
         margin-bottom: 0.45rem;
         max-width: 20ch;
         font-size: clamp(2rem, 4.1vw, 3.35rem);
         color: var(--blue-deep);
         line-height: 0.98;
         font-weight: 800;
         letter-spacing: -0.04em;
       }
       .hero-copy p {
         position: relative;
         max-width: 48rem;
         margin-bottom: 0;
         color: var(--muted);
         font-size: 1rem;
         line-height: 1.6;
       }
       .metric-card {
         min-height: 168px;
         border-radius: 24px;
         background: linear-gradient(180deg, rgba(255,255,255,0.98) 0%, rgba(248,251,255,0.95) 100%);
         border: 1px solid rgba(46, 93, 140, 0.10);
         box-shadow: 0 16px 38px rgba(31, 95, 153, 0.07);
       }
       .metric-card .card-body {
         display: flex;
         flex-direction: column;
         justify-content: space-between;
         gap: 0.55rem;
         padding: 1.15rem 1.2rem 1rem 1.2rem;
         min-height: 100%;
       }
       .metric-label {
         font-size: 0.85rem;
         text-transform: uppercase;
         letter-spacing: 0.06em;
         color: var(--muted);
       }
       .metric-value {
         display: block;
         margin: 0.15rem 0;
         font-size: clamp(2rem, 3vw, 2.75rem);
         font-weight: 800;
         background: linear-gradient(135deg, var(--viridis-deep), var(--viridis-bright));
         -webkit-background-clip: text;
         background-clip: text;
         color: transparent;
         line-height: 1.02;
         white-space: nowrap;
         overflow: visible;
       }
       .metric-note {
         color: var(--muted);
         margin-top: 0.3rem;
         margin-bottom: 0;
       }
       .bslib-page-fill {
         overflow-y: auto;
       }
       .bslib-sidebar-layout > .main {
         gap: 1.2rem;
       }
       .sidebar {
         background: linear-gradient(180deg, #ffffff 0%, #f8fbff 100%);
         border-right: 1px solid var(--border-soft);
         max-height: calc(100vh - 7rem);
         overflow-y: auto;
         scrollbar-width: thin;
       }
       .sidebar::-webkit-scrollbar,
       .dataTables_scrollBody::-webkit-scrollbar {
         width: 10px;
         height: 10px;
       }
       .sidebar::-webkit-scrollbar-thumb,
       .dataTables_scrollBody::-webkit-scrollbar-thumb {
         background: rgba(45, 112, 142, 0.35);
         border-radius: 999px;
       }
       .sidebar .form-group:last-child {
         margin-bottom: 1rem;
       }
       .card {
         border-radius: 22px;
         border: 1px solid var(--border-soft);
         box-shadow: none;
         overflow: hidden;
       }
       .card-header,
       .card-footer {
         background: linear-gradient(180deg, #ffffff 0%, #f7fbff 100%);
         border-color: var(--border-soft);
         color: var(--blue-deep);
         font-weight: 600;
       }
       .form-control,
       .selectize-input,
       .selectize-dropdown,
       .irs--shiny .irs-bar,
       .irs--shiny .irs-single {
         border-color: rgba(31, 95, 153, 0.18);
       }
       .form-control,
       .selectize-input {
         border-radius: 14px;
         box-shadow: none;
         min-height: 44px;
       }
       .selectize-dropdown {
         border-radius: 16px;
         box-shadow: none;
         overflow: hidden;
       }
       .selectize-dropdown .active {
         background: rgba(59, 130, 196, 0.14);
         color: var(--blue-deep);
       }
       .selectize-input.focus,
       .form-control:focus {
         border-color: rgba(31, 95, 153, 0.45);
         box-shadow: 0 0 0 3px rgba(59, 130, 196, 0.12);
       }
       .btn-default,
       .btn-primary {
         border-radius: 999px;
       }
       .table {
         margin-bottom: 0;
       }
       .dataTables_wrapper {
         padding-bottom: 0.2rem;
       }
       .dataTables_scroll {
         border: 1px solid rgba(46, 93, 140, 0.10);
         border-radius: 18px;
         overflow: hidden;
       }
       .dataTables_scrollHead,
       .dataTables_scrollBody {
         background: #ffffff;
       }
       .dataTables_wrapper table.dataTable {
         table-layout: fixed;
       }
       .dataTables_scrollBody {
         min-height: 18rem;
       }
       .dataTables_filter,
       .dataTables_length,
       .dataTables_info,
       .dataTables_paginate {
         padding-top: 0.55rem;
       }
       .dataTables_wrapper .dataTables_scrollHead thead tr:nth-child(2),
       .dataTables_wrapper .DTFC_Cloned thead tr:nth-child(2) {
         display: none;
       }
       table.dataTable thead .sorting::before,
       table.dataTable thead .sorting::after,
       table.dataTable thead .sorting_asc::before,
       table.dataTable thead .sorting_asc::after,
       table.dataTable thead .sorting_desc::before,
       table.dataTable thead .sorting_desc::after {
         display: none !important;
         content: none !important;
         opacity: 0 !important;
       }
       .dataTables_wrapper .dataTables_scrollHead table.dataTable thead th::before,
       .dataTables_wrapper .dataTables_scrollHead table.dataTable thead th::after,
       .dataTables_wrapper .DTFC_Cloned table.dataTable thead th::before,
       .dataTables_wrapper .DTFC_Cloned table.dataTable thead th::after {
         display: none !important;
         content: none !important;
         opacity: 0 !important;
       }
       .table > thead > tr > th {
         padding: 1rem 1rem 0.9rem 1rem;
         font-size: 0.95rem;
         font-weight: 800;
         letter-spacing: -0.01em;
         text-align: center;
         vertical-align: middle;
         white-space: normal !important;
         word-break: normal;
         overflow-wrap: anywhere;
         line-height: 1.35;
         min-height: 5.25rem;
         color: var(--blue-deep);
         border-bottom: 1px solid var(--border-soft);
         background: linear-gradient(180deg, rgba(248,251,255,0.98), rgba(238,246,255,0.92));
       }
       .dataTables_wrapper .dataTable thead th {
         background-image: none !important;
       }
       .dataTables_wrapper .dataTables_scrollHeadInner,
       .dataTables_wrapper .dataTables_scrollHeadInner table {
         width: 100% !important;
       }
       .dataTables_wrapper .dataTables_scrollHead table.dataTable thead th,
       .dataTables_wrapper .dataTables_scrollBody table.dataTable thead th {
         text-align: center;
         vertical-align: middle;
         white-space: normal !important;
         line-height: 1.35;
       }
       .table > tbody > tr > td {
         text-align: center;
         vertical-align: middle;
         border-top: 1px solid rgba(46, 93, 140, 0.08);
       }
       #survival_counts table th,
       #survival_counts table td {
         text-align: center;
         vertical-align: middle;
       }
       @media (max-width: 767px) {
         .hero-copy h1 {
           font-size: 1.7rem;
           max-width: 12ch;
         }
         .hero-banner {
           padding: 1.2rem 0 1rem 0;
         }
         .metric-card {
           min-height: 150px;
         }
         .sidebar {
           max-height: none;
         }
       }"
      ))
    ),
    tags$div(
      class = "hero-banner",
      tags$div(
        class = "hero-copy",
        tags$span("TCGA-KIRC Explorer"),
        tags$h1("Resumen clínico, genes y supervivencia"),
        tags$p(
          "Consulta de forma rápida los análisis principales de la cohorte TCGA-KIRC."
        )
      )
    )
  ),
  nav_panel(
    "Resumen",
    layout_column_wrap(
      width = 1/4,
      card(
        class = "metric-card",
        card_body(
          tags$div(class = "metric-label", "Pacientes con datos clínicos"),
          tags$div(class = "metric-value", format_metric(n_distinct(clinical_data$patient_id))),
          tags$p(class = "metric-note", "Cohorte clínica preparada para el análisis")
        )
      ),
      card(
        class = "metric-card",
        card_body(
          tags$div(class = "metric-label", "Muestras transcriptómicas"),
          tags$div(class = "metric-value", format_metric(nrow(sample_data))),
          tags$p(class = "metric-note", "Tumor primario y tejido normal")
        )
      ),
      card(
        class = "metric-card",
        card_body(
          tags$div(class = "metric-label", "Genes disponibles"),
          tags$div(class = "metric-value", format_metric(nrow(expression_matrix))),
          tags$p(class = "metric-note", paste("Assay utilizado:", expression_bundle$assay_name))
        )
      ),
      card(
        class = "metric-card",
        card_body(
          tags$div(class = "metric-label", "Tumores primarios"),
          tags$div(class = "metric-value", sum(sample_data$tissue_group == "Tumor primario")),
          tags$p(class = "metric-note", paste("Tejido normal:", sum(sample_data$tissue_group == "Tejido normal")))
        )
      )
    ),
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        selectInput(
          inputId = "overview_variable",
          label = "Variable clínica a resumir",
          choices = c(
            "Estadio AJCC" = "ajcc_pathologic_stage",
            "Grado tumoral" = "tumor_grade",
            "Sexo" = "sex",
            "Estado vital" = "vital_status"
          ),
          selected = "ajcc_pathologic_stage"
        )
      ),
      card(
        full_screen = TRUE,
        card_header("Distribución clínica"),
        plotOutput("overview_plot", height = 420)
      ),
      card(
        full_screen = TRUE,
        card_header("Vista previa de variables clínicas"),
        DTOutput("overview_table")
      )
    )
  ),
  nav_panel(
    "Genes",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        selectizeInput(
          inputId = "gene_symbol",
          label = "Gen de interés",
          choices = NULL,
          selected = default_gene,
          options = list(placeholder = "Escribe un símbolo génico", maxOptions = 150)
        ),
        selectInput(
          inputId = "gene_grouping",
          label = "Agrupar expresión por",
          choices = c(
            "Tipo de tejido" = "tissue_group",
            "Estadio AJCC" = "ajcc_pathologic_stage",
            "Grado tumoral" = "tumor_grade",
            "Sexo" = "sex",
            "Estado vital" = "vital_status"
          ),
          selected = "tissue_group"
        ),
        checkboxInput(
          inputId = "only_tumor_samples",
          label = "Mostrar solo muestras tumorales",
          value = FALSE
        )
      ),
      card(
        full_screen = TRUE,
        card_header(textOutput("gene_plot_title")),
        plotOutput("gene_plot", height = 500)
      ),
      card(
        full_screen = TRUE,
        card_header("Resumen por grupo"),
        tableOutput("gene_summary")
      )
    )
  ),
  nav_panel(
    "Supervivencia",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        selectInput(
          inputId = "survival_view",
          label = "Curva de supervivencia",
          choices = c(
            "Supervivencia global" = "global",
            "Estadio AJCC" = "ajcc",
            "Grado tumoral" = "grade",
            "Edad al diagnóstico" = "age",
            "Sexo biológico" = "sex",
            "Expresión del gen seleccionado" = "gene"
          ),
          selected = "global"
        ),
        conditionalPanel(
          condition = "input.survival_view === 'gene'",
          selectizeInput(
            inputId = "survival_gene",
            label = "Gen para el punto de corte mediano",
            choices = NULL,
            selected = default_gene,
            options = list(placeholder = "Escribe un símbolo génico", maxOptions = 150)
          )
        )
      ),
      card(
        full_screen = TRUE,
        card_header("Curvas de Kaplan-Meier"),
        plotOutput("survival_plot", height = 480),
        card_footer(textOutput("survival_caption"))
      ),
      card(
        full_screen = TRUE,
        card_header("Tamaño de los grupos analizados"),
        tableOutput("survival_counts")
      )
    )
  )
)

server <- function(input, output, session) {
  updateSelectizeInput(
    session = session,
    inputId = "gene_symbol",
    choices = gene_choices,
    selected = default_gene,
    server = TRUE
  )

  updateSelectizeInput(
    session = session,
    inputId = "survival_gene",
    choices = gene_choices,
    selected = default_gene,
    server = TRUE
  )

  overview_labels <- c(
    tissue_group = "Tipo de tejido",
    ajcc_pathologic_stage = "Estadio AJCC",
    tumor_grade = "Grado tumoral",
    age_group = "Edad al diagnóstico",
    sex = "Sexo",
    vital_status = "Estado vital"
  )

  gene_data <- reactive({
    req(input$gene_symbol)
    req(input$gene_symbol %in% rownames(expression_matrix))

    gene_df <- analysis_data %>%
      mutate(expression_value = as.numeric(expression_matrix[input$gene_symbol, sample_id]))

    if (isTRUE(input$only_tumor_samples)) {
      gene_df <- gene_df %>% filter(tissue_group == "Tumor primario")
    }

    gene_df %>%
      filter(!is.na(.data[[input$gene_grouping]]), !is.na(expression_value))
  })

  survival_context <- reactive({
    base_survival <- clinical_data %>%
      filter(!is.na(survival_time_months), survival_time_months > 0)

    view <- input$survival_view %||% "global"

    if (identical(view, "global")) {
      return(list(
        data = base_survival,
        group_variable = NULL,
        subtitle = "Cohorte clínica completa",
        legend_title = NULL,
        strata_levels = "Cohorte completa"
      ))
    }

    if (identical(view, "ajcc")) {
      surv_df <- base_survival %>%
        filter(ajcc_pathologic_stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV")) %>%
        mutate(
          survival_group = factor(
            dplyr::recode(
              ajcc_pathologic_stage,
              "Stage I" = "Estadio I",
              "Stage II" = "Estadio II",
              "Stage III" = "Estadio III",
              "Stage IV" = "Estadio IV"
            ),
            levels = c("Estadio I", "Estadio II", "Estadio III", "Estadio IV")
          )
        )

      return(list(
        data = surv_df,
        group_variable = "survival_group",
        subtitle = "Estratificación por estadio patológico AJCC",
        legend_title = "Estadio AJCC",
        strata_levels = levels(surv_df$survival_group)
      ))
    }

    if (identical(view, "grade")) {
      surv_df <- base_survival %>%
        filter(tumor_grade %in% c("G1", "G2", "G3", "G4")) %>%
        mutate(survival_group = factor(tumor_grade, levels = c("G1", "G2", "G3", "G4")))

      return(list(
        data = surv_df,
        group_variable = "survival_group",
        subtitle = "Estratificación por grado tumoral",
        legend_title = "Grado tumoral",
        strata_levels = levels(surv_df$survival_group)
      ))
    }

    if (identical(view, "age")) {
      surv_df <- base_survival %>%
        mutate(
          survival_group = case_when(
            age_at_diagnosis < 50 ~ "<50 años",
            age_at_diagnosis <= 70 ~ "50-70 años",
            age_at_diagnosis > 70 ~ ">70 años",
            TRUE ~ NA_character_
          ),
          survival_group = factor(survival_group, levels = c("<50 años", "50-70 años", ">70 años"))
        ) %>%
        filter(!is.na(survival_group))

      return(list(
        data = surv_df,
        group_variable = "survival_group",
        subtitle = "Estratificación por grupos de edad al diagnóstico",
        legend_title = "Edad",
        strata_levels = levels(surv_df$survival_group)
      ))
    }

    if (identical(view, "sex")) {
      surv_df <- base_survival %>%
        filter(sex %in% c("male", "female")) %>%
        mutate(
          survival_group = factor(
            dplyr::recode(sex, "female" = "Mujer", "male" = "Hombre"),
            levels = c("Mujer", "Hombre")
          )
        )

      return(list(
        data = surv_df,
        group_variable = "survival_group",
        subtitle = "Estratificación por sexo biológico",
        legend_title = "Sexo",
        strata_levels = levels(surv_df$survival_group)
      ))
    }

    surv_df <- build_survival_dataset(
      clinical_data = clinical_data,
      sample_data = sample_data,
      expression_matrix = expression_matrix,
      gene_symbol = input$survival_gene
    ) %>%
      filter(!is.na(survival_time_months), survival_time_months > 0) %>%
      mutate(survival_group = factor(gene_group, levels = c("Alta expresión", "Baja expresión")))

    list(
      data = surv_df,
      group_variable = "survival_group",
      subtitle = paste("Estratificación por expresión mediana de", input$survival_gene),
      legend_title = "Expresión",
      strata_levels = levels(surv_df$survival_group)
    )
  })

  output$overview_plot <- renderPlot({
    overview_df <- clinical_data %>%
      transmute(group = as.character(.data[[input$overview_variable]])) %>%
      filter(!is.na(group), nzchar(group)) %>%
      group_by(group) %>%
      summarise(n = dplyr::n(), .groups = "drop")

    ggplot(overview_df, aes(x = reorder(group, n), y = n, fill = group)) +
      geom_col(width = 0.75, show.legend = FALSE) +
      coord_flip() +
      scale_fill_viridis_d(option = "D", begin = 0.34, end = 0.62) +
      labs(
        title = paste("Distribución de", tolower(get_variable_label(input$overview_variable, overview_labels))),
        x = NULL,
        y = "Número de pacientes"
      ) +
      theme(
        axis.text.y = element_text(size = 11),
        plot.margin = margin(8, 18, 8, 6)
      )
  })

  output$overview_table <- renderDT({
    clinical_data %>%
      dplyr::select(
        patient_id,
        age_at_diagnosis,
        ajcc_pathologic_stage,
        tumor_grade,
        sex,
        vital_status,
        survival_time_months
      ) %>%
      dplyr::rename(
        `ID del paciente` = patient_id,
        `Edad al diagnóstico` = age_at_diagnosis,
        `Estadio AJCC` = ajcc_pathologic_stage,
        `Grado tumoral` = tumor_grade,
        Sexo = sex,
        `Estado vital` = vital_status,
        `Supervivencia (meses)` = survival_time_months
      ) %>%
      datatable(
        rownames = FALSE,
        filter = "none",
        class = "stripe hover compact",
        options = list(
          pageLength = 12,
          lengthMenu = c(12, 25, 50, 100),
          scrollX = TRUE,
          scrollY = "420px",
          scrollCollapse = TRUE,
          fixedHeader = TRUE,
          autoWidth = FALSE,
          columnDefs = list(
            list(width = "13%", targets = 0),
            list(width = "12%", targets = 1),
            list(width = "12%", targets = 2),
            list(width = "12%", targets = 3),
            list(width = "8%", targets = 4),
            list(width = "11%", targets = 5),
            list(width = "14%", targets = 6)
          ),
          searching = TRUE
        ),
        callback = JS(
          "table.on('init.dt draw.dt', function() {",
          "  var container = $(table.table().container());",
          "  container.find('.dataTables_scrollHead thead tr:gt(0)').remove();",
          "  container.find('thead th').css('background-image', 'none');",
          "});"
        )
      ) %>%
      formatRound(columns = c("Edad al diagnóstico", "Supervivencia (meses)"), digits = 2, dec.mark = ",")
  })

  output$gene_plot_title <- renderText({
    paste("Expresión de", input$gene_symbol)
  })

  output$gene_plot <- renderPlot({
    gene_df <- gene_data() %>%
      mutate(group = factor(.data[[input$gene_grouping]]))

    group_label <- tolower(get_variable_label(input$gene_grouping, overview_labels))

    ggplot(gene_df, aes(x = group, y = log2(expression_value + 1))) +
      geom_boxplot(aes(fill = group), alpha = 0.94, outlier.shape = NA, width = 0.48, show.legend = FALSE, colour = "#355f8d", linewidth = 0.55) +
      geom_point(
        position = position_jitter(width = 0.12, height = 0),
        alpha = 0.24,
        size = 1.25,
        colour = "#2d708e",
        show.legend = FALSE
      ) +
      scale_fill_viridis_d(option = "D", begin = 0.12, end = 0.88) +
      labs(
        title = paste("Expresión de", input$gene_symbol, "según", group_label),
        subtitle = paste("Valores representados como log2(", expression_bundle$assay_name, "+ 1)", sep = ""),
        x = NULL,
        y = "Expresión transformada"
      ) +
      theme(axis.text.x = element_text(angle = 18, hjust = 1))
  })

  output$gene_summary <- renderTable({
    gene_data() %>%
      mutate(group = .data[[input$gene_grouping]]) %>%
      group_by(group) %>%
      summarise(
        muestras = n(),
        media = round(mean(log2(expression_value + 1), na.rm = TRUE), 2),
        mediana = round(median(log2(expression_value + 1), na.rm = TRUE), 2),
        .groups = "drop"
      )
  }, striped = TRUE, bordered = FALSE, spacing = "s")

  output$survival_plot <- renderPlot({
    ctx <- survival_context()
    surv_df <- ctx$data
    req(nrow(surv_df) > 1)

    fit <- if (is.null(ctx$group_variable)) {
      survfit(
        Surv(survival_time_months, death_event) ~ 1,
        data = surv_df
      )
    } else {
      survfit(
        as.formula(paste("Surv(survival_time_months, death_event) ~", ctx$group_variable)),
        data = surv_df
      )
    }

    km_data <- make_km_data(fit)
    req(nrow(km_data) > 0)

    if (!is.null(ctx$strata_levels)) {
      km_data <- km_data %>%
        mutate(strata = factor(strata, levels = ctx$strata_levels))
    }

    survival_plot <- ggplot(km_data, aes(x = time, y = surv)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        title = "Curvas de Kaplan-Meier",
        subtitle = ctx$subtitle,
        x = "Tiempo de seguimiento (meses)",
        y = "Supervivencia estimada"
      )

    if (is.null(ctx$group_variable)) {
      survival_plot +
        geom_step(colour = "#2D708E", linewidth = 1.15)
    } else {
      survival_plot +
        geom_step(aes(colour = strata), linewidth = 1.1) +
        scale_colour_viridis_d(option = "D", begin = 0.30, end = 0.72, drop = FALSE) +
        labs(colour = ctx$legend_title)
    }
  })

  output$survival_caption <- renderText({
    ctx <- survival_context()
    surv_df <- ctx$data
    req(nrow(surv_df) > 1)

    if (is.null(ctx$group_variable)) {
      fit <- survfit(Surv(survival_time_months, death_event) ~ 1, data = surv_df)
      fit_table <- summary(fit)$table
      median_value <- unname(fit_table["median"])

      if (is.na(median_value)) {
        return("Supervivencia global de la cohorte clínica")
      }

      return(paste("Mediana de supervivencia estimada:", round(median_value, 1), "meses"))
    }

    diff_object <- survdiff(
      as.formula(paste("Surv(survival_time_months, death_event) ~", ctx$group_variable)),
      data = surv_df
    )

    p_value <- 1 - pchisq(diff_object$chisq, df = max(length(diff_object$n) - 1, 1))
    paste("Prueba de log-rank: p =", format.pval(p_value, digits = 3, eps = 0.001))
  })

  output$survival_counts <- renderTable({
    ctx <- survival_context()
    surv_df <- ctx$data
    req(nrow(surv_df) > 1)

    if (is.null(ctx$group_variable)) {
      return(
        tibble::tibble(
          grupo = "Cohorte completa",
          pacientes = nrow(surv_df),
          eventos = sum(surv_df$death_event, na.rm = TRUE)
        )
      )
    }

    group_counts <- surv_df %>%
      mutate(grupo = as.character(.data[[ctx$group_variable]])) %>%
      filter(!is.na(grupo), nzchar(grupo)) %>%
      group_by(grupo) %>%
      summarise(
        pacientes = dplyr::n(),
        eventos = sum(death_event, na.rm = TRUE),
        .groups = "drop"
      )

    if (!is.null(ctx$strata_levels)) {
      group_counts <- group_counts %>%
        mutate(grupo = factor(grupo, levels = ctx$strata_levels)) %>%
        arrange(grupo) %>%
        mutate(grupo = as.character(grupo))
    }

    group_counts
  }, striped = TRUE, bordered = FALSE, spacing = "s", align = "ccc")
}

app <- shinyApp(ui = ui, server = server)

if (sys.nframe() == 0) {
  runApp(app)
}

app