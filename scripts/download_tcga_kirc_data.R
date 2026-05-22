suppressPackageStartupMessages({
  library(TCGAbiolinks)
})

project_dir <- getwd()

if (basename(project_dir) == "scripts") {
  project_dir <- dirname(project_dir)
}

data_dir <- file.path(project_dir, "data")
transcriptome_dir <- file.path(data_dir, "transcriptome")
clinical_dir <- file.path(data_dir, "clinical")
clinical_file <- file.path(clinical_dir, "clinical_data_TCGA_KIRC.tsv")

dir.create(transcriptome_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clinical_dir, recursive = TRUE, showWarnings = FALSE)

transcriptome_files <- list.files(
  transcriptome_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(transcriptome_files) == 0) {
  query_exp <- GDCquery(
    project = "TCGA-KIRC",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  GDCdownload(
    query = query_exp,
    method = "api",
    files.per.chunk = 20,
    directory = transcriptome_dir
  )
} else {
  message("Los archivos transcriptómicos ya están disponibles.")
}

if (!file.exists(clinical_file)) {
  clinical_data <- GDCquery_clinic(
    project = "TCGA-KIRC",
    type = "clinical"
  )

  write.table(
    clinical_data,
    file = clinical_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
} else {
  message("El archivo clínico ya está disponible.")
}