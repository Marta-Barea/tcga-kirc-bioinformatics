suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

input_path <- "data/prepared/tcga_kirc_star_counts_se.rds"
output_path <- "data/prepared/tcga_kirc_star_counts_se_small.rds"

candidate_genes <- c(
  "CA9", "NDUFA4L2", "EGLN3", "HILPDA", "SCARB1", "STC2",
  "COL23A1", "CDKN2A", "GABRD", "PVT1", "TTC21B-AS1", "LINC00887"
)

message("Leyendo archivo grande...")
se <- readRDS(input_path)

gene_names <- rowData(se)$gene_name
if (is.null(gene_names)) {
  gene_names <- rownames(se)
}

keep <- gene_names %in% candidate_genes | rownames(se) %in% candidate_genes

message("Genes encontrados: ", sum(keep))

if (sum(keep) == 0) {
  stop("No se encontró ningún gen candidato en el objeto.")
}

se_small <- se[keep, ]

message("Guardando archivo reducido en: ", output_path)
saveRDS(se_small, output_path, compress = "xz")

message("Tamaño original MB:")
print(file.info(input_path)$size / 1024^2)

message("Tamaño reducido MB:")
print(file.info(output_path)$size / 1024^2)

message("Listo.")