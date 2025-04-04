#######
# Script: Visualización de resultados de expresión diferencial (Volcano Plot y Heatmaps)
# Author: Paola Godoy, Yael Hernández y Ariadna Badía
# Date: 13/03/2025
# Description: Este script genera visualizaciones de los resultados de expresión diferencial,
# incluyendo un Volcano Plot para mostrar genes regulados positivamente y negativamente,
# y Heatmaps para visualizar los genes más significativos y sus cambios en expresión.
# Arguments:
#   - Input: Resultados de expresión diferencial (DE_ATF3_vs_vehiculo.csv) (DE_Aart6_vs_ATF3.csv),
#            objetos dds y vsd (dds_Times_vs_control.RData y vst_Times_vs_control.RData).
#   - Output: Volcano Plot y Heatmaps (archivos PNG).
#######

######
#Heatmap
#Propósito: Mostrar patrones de expresión en todas las condiciones (Vehículo, ATF3, Aart6), sin comparaciones directas.
#Base de datos: Utiliza los valores normalizados (rlog o vst) de todas las muestras, no los resultados de expresión diferencial (DEG).

# --- Load packages ----------
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(org.Hs.eg.db)

# --- Load data -----
# Cargar variable "dds", proveniente del script "DEG_analysis.R"
figdir = "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/figures/"

# Cargar variable "dds", proveniente del script "DEG_analysis.R"
load("/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/dds_Times_vs_vehiculo.RData")

# Cargar variable "vsdata", proveniente del script "DEG_analysis.R"
load("/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/vst_Times_vs_vehiculo.RData")

# 1. Cargar datos normalizados (rlog)
rld <- rlog(dds, blind = FALSE)
mat <- assay(rld)

# 2. Definir genes de interés (extraídos del script original)
genes <- c(
  "IFNG","IFNGR1", "IFNGR2", "TNF", "HLA-A", "HLA-B", "HLA-C",
  "HLA-F", "HLA-J", "HLA-DPB2", "HLA-DRA", "HLA-DQA2","IL11", 
  "IL4","SLCO3A1","AIF1L", "CRABP2", "IRF8", "IRF9","TAP1", "STAT1", "STAT2", 
  "STAT3","STAT4","STAT5A", "STAT5B","STAT6","MAF", "REL", "RELA",
  "RELB","NFKB1","NFKB2"
)

# 3. Convertir símbolos a IDs de Ensembl (sin versión)
ensembl_ids <- mapIds(
  org.Hs.eg.db,
  keys = genes,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)
ensembl_ids <- na.omit(ensembl_ids)

# 4. Filtrar matriz y ajustar nombres
rownames(mat) <- sub("\\..*", "", rownames(mat)) # Eliminar versiones
mat_filtrado <- mat[rownames(mat) %in% ensembl_ids, ]
rownames(mat_filtrado) <- names(ensembl_ids)[match(rownames(mat_filtrado), ensembl_ids)]

# 5. Ordenar las muestras como en el artículo: Vehículo -> ATF3 -> Aart6
metadata_ordenado <- metadata[order(metadata$type, decreasing = FALSE), ]
mat_filtrado <- mat_filtrado[, metadata_ordenado$sample_id]

# 6. Anotaciones y colores (extraídos del script original)
annotation_col <- data.frame(
  Condition = factor(
    metadata_ordenado$type,
    levels = c("vehiculo", "ATF3", "Aart6")
  ),
  row.names = metadata_ordenado$sample_id
)
ann_colors <- list(
  Condition = c(
    vehiculo = "purple",
    ATF3 = "#F5A",
    Aart6 = "#7ED321"      # Verde
  )
)

# 7. Parámetros estéticos exactos del artículo
png(
  file = paste0(figdir, "Heatmap.png"),
  width = 1800,  # Ancho en píxeles (alta resolución)
  height = 1500,
  res = 300      # 300 DPI para publicación
)
pheatmap(
  mat_filtrado,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), # Mismo esquema de color
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # Muestras ordenadas manualmente
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_colnames = TRUE,  # Ocultar nombres de muestras
  fontsize_row = 9,
  fontsize_col = 12,
  border_color = NA,
  gaps_col = cumsum(table(metadata_ordenado$type)), # Separadores entre condiciones
  main = ""
)
dev.off()


#####
#Volcano plot
# --- Código para Volcano Plots con genes destacados ----
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(patchwork)

# 1. Función para generar Volcano Plots con genes clave
generate_volcano <- function(res_path, contrast_name, key_genes) {
  # Cargar resultados DEG
  res <- read.csv(res_path, row.names = 1)
  
  # --- Limpiar IDs de Ensembl ----
  rownames(res) <- sub("\\..*", "", rownames(res))
  
  # Convertir IDs de Ensembl a símbolos génicos
  ensembl_ids <- rownames(res)
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Crear dataframe con símbolos
  df <- as.data.frame(res) %>%
    mutate(
      Gene = gene_symbols,
      Expression = case_when(
        log2FoldChange >= 2 & padj < 0.05 ~ "Up-regulated",
        log2FoldChange <= -2 & padj < 0.05 ~ "Down-regulated",
        TRUE ~ "Unchanged"
      ),
      KeyGene = ifelse(Gene %in% key_genes, "Yes", "No")
    ) %>%
    filter(!is.na(Gene))  
  
  # Filtrar genes clave significativos
  df_key <- df %>% filter(KeyGene == "Yes", padj < 0.05, abs(log2FoldChange) >= 2)
  
  # Generar plot
  p <- ggplot(df, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(color = Expression), size = 0.7, alpha = 0.6) +
    geom_text_repel(
      data = df_key,
      aes(label = Gene),
      color = "black",
      size = 3,
      box.padding = 0.5,
      max.overlaps = 50
    ) +
    scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
    labs(
      title = contrast_name,
      x = expression("log"[2]*"FC"),
      y = expression("-log"[10]*"p-adj")
    ) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    theme_minimal()
  
  return(p)
}

# 2. Definir genes clave del artículo
key_genes <- c(
  "IFNG", "IFNGR1", "IFNGR2", "TNF", "HLA-A", "HLA-B", "HLA-C",
  "HLA-F", "HLA-DRA", "IL11", "SLCO3A1", "AIF1L", "CRABP2", 
  "IRF8", "IRF9", "TAP1", "STAT1", "STAT2", "STAT3", "STAT4", 
  "STAT5A", "STAT5B", "STAT6", "MAF", "REL", "RELA", "RELB", 
  "NFKB1", "NFKB2"
)

# 3. Generar Volcano Plots para cada contraste

p1 = generate_volcano(
  "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/DE_ATF3_vs_vehiculo.csv",
  "ATF3 vs Vehículo",
  key_genes
)

p2 = generate_volcano(
  "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/DE_Aart6_vs_ATF3.csv",
  "ATF3 vs Aart6",
  key_genes
)

p3 = generate_volcano(
  "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/DE_Aart6_vs_vehiculo.csv",
  "Aart6 vs Vehículo",
  key_genes
)

# Combinar en una sola imagen (1 fila x 3 columnas)
combined_plot <- p1 + p2 + p3 + 
  plot_layout(ncol = 3) + 
  plot_annotation(tag_levels = "A")  # Etiquetas A, B, C

# Guardar la imagen combinada
ggsave(
  paste0(figdir, "Combined_Volcano.png"),
  combined_plot,
  width = 18,  # Ancho total ajustado para 3 columnas
  height = 6,   # Altura constante
  dpi = 300
)





