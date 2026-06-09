[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R](https://img.shields.io/badge/R-Bioinformatics-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Shiny App](https://img.shields.io/badge/Shiny-Live_App-1f5f99?logo=rstudioide&logoColor=white)](https://marta-barea-sepul.shinyapps.io/tcga-kirc-explorer/)
[![Report PDF](https://img.shields.io/badge/Report-PDF-1b5e20)](reports/tcga_kirc_bioinformatics.pdf)
[![TCGAbiolinks](https://img.shields.io/badge/TCGAbiolinks-Bioconductor-0f7b0f)](https://bioconductor.org/packages/TCGAbiolinks/)

# TCGA-KIRC Bioinformatics

## Overview

This repository implements a reproducible R-based workflow for the analysis of TCGA-KIRC (The Cancer Genome Atlas Kidney Renal Clear Cell Carcinoma) data. The project integrates clinical and transcriptomic data, performs exploratory and inferential analyses, generates a formal report in PDF format, and exposes a lightweight Shiny application for interactive inspection of selected candidate genes.

The analytical scope covers:

- acquisition of TCGA-KIRC clinical and RNA-seq STAR counts data from the Genomic Data Commons;
- preparation of serialized R objects for downstream analyses;
- exploratory clinical profiling and survival analysis;
- transcriptomic preprocessing, normalization, variance-stabilizing transformation, PCA, and differential expression analysis;
- integration of candidate gene expression with clinical survival outcomes;
- interactive visualization through a Shiny dashboard.

## Repository Components

### Data ingestion and preparation

The project downloads TCGA-KIRC resources directly from GDC using TCGAbiolinks and stores prepared artifacts as RDS files under the data/prepared directory.

- scripts/download_tcga_kirc_data.R: downloads transcriptomic and clinical datasets and creates the prepared R objects.
- scripts/create_small_star_counts.R: generates a reduced SummarizedExperiment containing a focused candidate-gene panel for lightweight interactive use.

### Analytical report

The main report is authored in R Markdown and rendered to PDF.

- reports/tcga_kirc_bioinformatics.Rmd: end-to-end analytical document.
- reports/tcga_kirc_bioinformatics.pdf: rendered technical report.
- reports/references.bib and reports/preamble.tex: bibliography and LaTeX customization assets.

### Interactive application

The Shiny application provides exploratory access to the prepared clinical dataset and the reduced transcriptomic object used for candidate-gene inspection.

- app/app.R: standalone Shiny application.
- deployed app: https://marta-barea-sepul.shinyapps.io/tcga-kirc-explorer/

## Project Structure

```text
.
├── app/                         # Shiny application and deployment metadata
├── assets/                      # Static assets used by the project
├── cache/                       # Knitr/cache artifacts from report execution
├── data/
│   ├── clinical/                # Raw downloaded clinical tables
│   ├── prepared/                # Serialized R objects used by report and app
│   └── transcriptome/           # Downloaded TCGA transcriptomic files
├── reports/                     # R Markdown source, bibliography, LaTeX assets, PDF output
├── scripts/                     # Data acquisition and preprocessing scripts
├── LICENSE
├── MANIFEST.txt                 # GDC manifest used for transcriptomic downloads
└── README.md
```

## Data Sources

This workflow uses publicly available TCGA-KIRC resources from the Genomic Data Commons:

- clinical data obtained through the TCGAbiolinks clinical query interface;
- transcriptome profiling data from the Transcriptome Profiling category;
- gene expression quantification generated with the STAR - Counts workflow.

The manifest included in the repository documents the downloaded transcriptomic files associated with the analysis.

## Main Outputs

The repository produces the following main artifacts:

- prepared clinical object: data/prepared/clinical_data_TCGA_KIRC.rds;
- prepared transcriptomic object: data/prepared/tcga_kirc_star_counts_se.rds;
- reduced transcriptomic object for the app: data/prepared/tcga_kirc_star_counts_se_small.rds;
- rendered analytical report: reports/tcga_kirc_bioinformatics.pdf.

## Requirements

The project is implemented in R and relies on CRAN and Bioconductor packages. Core dependencies inferred from the report, scripts, and app include:

- TCGAbiolinks
- SummarizedExperiment
- DESeq2
- survival
- shiny
- bslib
- DT
- dplyr
- ggplot2
- magrittr
- kableExtra

A LaTeX engine is also required to build the PDF report, as the R Markdown configuration uses xelatex.

## Setup

Clone the repository and install the required packages in your R environment.

```bash
git clone git@github.com:Marta-Barea/tcga-kirc-bioinformatics.git
cd tcga-kirc-bioinformatics
```

Example package installation in R:

```r
install.packages(c(
	"shiny", "bslib", "DT", "dplyr", "ggplot2", "magrittr", "survival", "kableExtra"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
	install.packages("BiocManager")
}

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2"))
```

## Reproducible Workflow

### 1. Download and prepare TCGA-KIRC data

Run the acquisition script from the repository root:

```bash
Rscript scripts/download_tcga_kirc_data.R
```

This step:

- queries TCGA-KIRC transcriptomic STAR counts data from GDC;
- downloads the files into data/transcriptome;
- prepares a SummarizedExperiment object and stores it in data/prepared;
- downloads clinical annotations and stores them both as TSV and RDS.

### 2. Build the reduced transcriptomic object for the Shiny app

```bash
Rscript scripts/create_small_star_counts.R
```

This script extracts a compact candidate-gene panel from the full transcriptomic object to reduce application footprint during interactive use and deployment.

### 3. Render the technical report

```bash
Rscript -e "rmarkdown::render('reports/tcga_kirc_bioinformatics.Rmd', quiet = TRUE)"
```

The rendered PDF is written to reports/tcga_kirc_bioinformatics.pdf.

### 4. Launch the Shiny app locally

```bash
Rscript -e "shiny::runApp('app')"
```

## Analytical Scope

The report and the application cover the following technical tasks:

- clinical data curation and missingness handling;
- descriptive analysis of categorical and continuous clinical variables;
- Kaplan-Meier and Cox proportional hazards modeling;
- inspection of transcriptomic count matrices and gene identifier conversion;
- DESeq2-based normalization and differential expression analysis;
- variance-stabilizing transformation and PCA for transcriptomic structure assessment;
- focused evaluation of biologically relevant candidate genes, including CA9, NDUFA4L2, EGLN3, HILPDA, SCARB1, STC2, COL23A1, CDKN2A, GABRD, PVT1, TTC21B-AS1, and LINC00887.

## Notes on Reproducibility

- The data directory is intentionally ignored by Git because downloaded TCGA artifacts and prepared R objects are large.
- The cache directory contains intermediate knitr artifacts generated during report execution.
- The Shiny app expects the prepared clinical object and the reduced transcriptomic object to exist before launch.
- The report and app are designed to be executed from the repository root or from subdirectories that resolve back to the project root.

## Access

- Interactive app: https://marta-barea-sepul.shinyapps.io/tcga-kirc-explorer/
- PDF report: [reports/tcga_kirc_bioinformatics.pdf](reports/tcga_kirc_bioinformatics.pdf)

## Author

Marta Barea Sepulveda  
PhD  
Interuniversity Master's Degree in Bioinformatics and Biostatistics  
Universitat Oberta de Catalunya - Universitat de Barcelona

## License

This project is distributed under the GNU General Public License v3.0. See the LICENSE file for the full license text.
