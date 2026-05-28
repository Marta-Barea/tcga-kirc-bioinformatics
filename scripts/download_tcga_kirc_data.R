# Carga las librerias necesarias
suppressPackageStartupMessages({
  library(TCGAbiolinks)
})

# Función para eliminar archivos si existen
cleanup_files <- function(paths) {
  existing_paths <- paths[file.exists(paths)]

  if (length(existing_paths) > 0) {
    invisible(file.remove(existing_paths))
  }
}

# Establece los directorios y rutas de archivos
project_dir <- getwd()

if (basename(project_dir) == "scripts") {
  project_dir <- dirname(project_dir)
}

data_dir <- file.path(project_dir, "data")
transcriptome_dir <- file.path(data_dir, "transcriptome")
clinical_dir <- file.path(data_dir, "clinical")
prepared_dir <- file.path(data_dir, "prepared")
clinical_file <- file.path(clinical_dir, "clinical_data_TCGA_KIRC.tsv")
prepared_file <- file.path(prepared_dir, "tcga_kirc_star_counts_se.rds")
prepared_clinical_file <- file.path(prepared_dir, "clinical_data_TCGA_KIRC.rds")
auxiliary_rds_root <- file.path(project_dir, c("df.rds", "results.rds"))
auxiliary_rds_prepared <- file.path(prepared_dir, c("df.rds", "results.rds"))

# Crea los directorios necesarios si no existen
dir.create(transcriptome_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clinical_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)

# Define la consulta para los datos transcriptómicos
query_exp <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

transcriptome_files <- list.files(
  transcriptome_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# Descarga los datos transcriptómicos si no están disponibles
if (length(transcriptome_files) == 0) {
  GDCdownload(
    query = query_exp,
    method = "api",
    files.per.chunk = 20,
    directory = transcriptome_dir
  )
} else {
  message("Los archivos transcriptómicos ya están disponibles.")
}

# Prepara los datos transcriptómicos y guarda los objetos preparados
if (!file.exists(prepared_file)) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(prepared_dir)

  prepared_data <- GDCprepare(
    query = query_exp,
    directory = transcriptome_dir
  )

  cleanup_files(auxiliary_rds_prepared)
  cleanup_files(auxiliary_rds_root)

  saveRDS(prepared_data, file = prepared_file)
} else {
  message("El objeto transcriptómico preparado ya está disponible.")
}

# Descarga los datos clínicos si no están disponibles
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

# Prepara los datos clínicos y carda el objeto preparado
if (!file.exists(prepared_clinical_file)) {
  if (!exists("clinical_data")) {
    clinical_data <- read.delim(
      clinical_file,
      sep = "\t",
      check.names = FALSE
    )
  }

  saveRDS(clinical_data, file = prepared_clinical_file)
} else {
  message("El objeto clínico preparado ya está disponible.")
}

cleanup_files(auxiliary_rds_root)