#######
# Script: Análisis de enriquecimiento funcional y visualización de términos GO
# Author: Paola Godoy, Yael Hernández y Ariadna Badía
# Date: 15/03/2025
# Description: Este script realiza un análisis de enriquecimiento funcional utilizando la herramienta gprofiler2
# para genes regulados positivamente (up-regulated) y negativamente (down-regulated). Además, genera gráficos
# de Manhattan y barplots para visualizar los términos GO y otras categorías enriquecidas.
# Arguments:
#   - Input: Archivos CSV con resultados de expresión diferencial (log2FoldChange y padj).
#   - Output: Gráficos de Manhattan y barplots, así como archivos RData con los resultados del enriquecimiento.
#######


# --- Load packages ----------
library(gprofiler2)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(dplyr)

# --- Load data -----
# Cargar archivos
indir <- "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/"
outdir <- "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/"
figdir <- "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/figures/"

# Crear directorio de figuras si no existe
if (!dir.exists(figdir)) {
  dir.create(figdir, recursive = TRUE)
}

# Seleccionar todos los archivos CSV
files <- dir(indir, pattern = "\\.csv$")  # Selecciona todos los archivos CSV

# Verificar los archivos encontrados
print(files)  # Muestra los archivos CSV que se van a procesar

# ---- Bucle para procesar cada archivo CSV ----
for (file in files) {
  # Extraer el nombre del archivo sin extensión
  plot_name <- gsub("\\.csv$", "", file)  # Elimina la extensión .csv para usar en nombres de salida

  # Cargar archivo
  df <- read.csv(file = paste0(indir, file), row.names = 'X')
  head(df)

  # Agregar informacion sobre la expresion
  abslogFC <- 1  # Umbral más relajado para log2FoldChange
  df <- df %>%
    dplyr::mutate(Expression = case_when(log2FoldChange >= abslogFC & padj < 0.5 ~ "Up-regulated",
                                         log2FoldChange <= -(abslogFC) & padj < 0.5 ~ "Down-regulated",
                                         TRUE ~ "Unchanged"))

  # Obtener los nombres de los genes
  # > UP
  up_genes <- df %>% filter(Expression == 'Up-regulated') %>%
    arrange(padj, desc(abs(log2FoldChange)))
  # Extraer solo el nombre de los genes
  up_genes <- rownames(up_genes)

  # > Down
  down_genes <- df %>% filter(Expression == 'Down-regulated') %>%
    arrange(padj, desc(abs(log2FoldChange)))
  # Extraer solo el nombre de los genes
  down_genes <- rownames(down_genes)

  # --- Eliminar la parte de la versión de los identificadores de Ensembl ---
  up_genes <- gsub("\\..*", "", up_genes)
  down_genes <- gsub("\\..*", "", down_genes)

  # --- Imprimir los genes de entrada para verificar ---
  print("Up-regulated genes:")
  print(up_genes)

  print("Down-regulated genes:")
  print(down_genes)

  # --- Definir Category_colors ---
  Category_colors <- data.frame(
    category = c("GO:BP", "GO:CC", "GO:MF", "KEGG",
                 'REAC', 'TF', 'MIRNA', 'HPA', 'CORUM', 'HP', 'WP'),
    label = c('Biological Process', 'Cellular Component', 'Molecular Function',  "KEGG",
              'REAC', 'TF', 'MIRNA', 'HPA', 'CORUM', 'HP', 'WP'),
    colors =  c('#FF9900', '#109618','#DC3912', '#DD4477',
                '#3366CC','#5574A6', '#22AA99', '#6633CC', '#66AA00', '#990099', '#0099C6')
  )

  # ---- Analisis de terminos Go ----
  # Seleccionar bases de datos
  sources_db <- c("GO:BP", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP")

  # Realizar el enriquecimiento funcional
  multi_gp <- gost(list("Upregulated" = up_genes,
                        "Downregulated" = down_genes),
                   correction_method = "fdr", user_threshold = 0.5,  # Umbral
                   multi_query = F, ordered_query = T,
                   sources = sources_db,
                   evcodes = TRUE,
                   organism = 'hsapiens')  # Asegúrate de que el organismo sea correcto

  # Verificar si hay resultados
  if (!is.null(multi_gp$result)) {
    ## ----manhattan plot--------
    gostp1 <- gostplot(multi_gp, interactive = FALSE)
    # Guardar grafica
    ggsave(paste0(figdir, "ManhattanGO_", plot_name, ".png"),
           plot = gostp1, dpi = 300)

    ## ----Dataframe de todos los datos --------
    # Convertir a dataframe
    gost_query <- as.data.frame(multi_gp$result)

    # Extarer informacion en modo matriz de todos los resultados
    bar_data <- data.frame("term" = as.factor(gost_query$term_name), "condition" = gost_query$query,
                           "count" = gost_query$term_size, "p.adjust" = gost_query$p_value,
                           'category' = as.factor(gost_query$source), "go_id" = as.factor(gost_query$term_id),
                           'geneNames' = gost_query$intersection
    )

    ## ---- DOWN genes ----
    bar_data_down <- subset(bar_data, condition == 'Downregulated')

    # Verificar si hay datos en bar_data_down
    if (nrow(bar_data_down) > 0) {
      # Ordenar datos y seleccion por pvalue
      bar_data_down <- head(bar_data_down[order(bar_data_down$p.adjust),], 40) # order by pvalue
      bar_data_down_ordered <- bar_data_down[order(bar_data_down$p.adjust),] # order by pvalue
      bar_data_down_ordered <- bar_data_down_ordered[order(bar_data_down_ordered$category),] # order by category
      bar_data_down_ordered$p.val <- round(-log10(bar_data_down_ordered$p.adjust), 2)
      bar_data_down_ordered$num <- seq(1:nrow(bar_data_down_ordered)) # num category for plot

      # Guardar dataset
      save(bar_data_down_ordered, file = paste0(outdir, "DOWN_GO_", plot_name, ".RData"))

      # agregar colores para la grafica
      bar_data_down_ordered_mod <- left_join(bar_data_down_ordered, Category_colors, by= "category")

      ### ---- DOWN genes (barplot) ----
      # Generar la grafica
      g.down <- ggplot(bar_data_down_ordered_mod, aes(p.val, reorder(term, -num), fill = category)) +
        geom_bar(stat = "identity") +
        geom_text(
          aes(label = p.val),
          color = "black",
          hjust = 0,
          size = 2.2,
          position = position_dodge(0)
        ) +
        labs(x = "-log10(p-value)", y = NULL) +
        scale_fill_manual(
          name = 'Category',
          labels = unique(bar_data_down_ordered_mod$label),
          values = unique(bar_data_down_ordered_mod$colors)
        ) +
        theme(
          legend.position = "right",
          axis.title.y = element_blank(),
          strip.text.x = element_text(size = 11, face = "bold"),
          strip.background = element_blank()
        ) + theme_classic()

      # Guardar la figura
      ggsave(paste0(figdir, "barplotDOWN_GO_", plot_name, ".png"),
             plot = g.down + theme_classic(), dpi = 600, width = 10, height = 5)
    } else {
      message("No hay términos GO significativos para genes down-regulados en ", plot_name)
    }

    ## ---- UP genes ----
    bar_data_up <- subset(bar_data, condition == 'Upregulated')

    # Verificar si hay datos en bar_data_up
    if (nrow(bar_data_up) > 0) {
      # Ordenar datos y seleccion por pvalue
      bar_data_up <- head(bar_data_up[order(bar_data_up$p.adjust),], 40) # order by pvalue
      bar_data_up_ordered <- bar_data_up[order(bar_data_up$p.adjust),] # order by pvalue
      bar_data_up_ordered <- bar_data_up_ordered[order(bar_data_up_ordered$category),] # order by category
      bar_data_up_ordered$p.val <- round(-log10(bar_data_up_ordered$p.adjust), 2)
      bar_data_up_ordered$num <- seq(1:nrow(bar_data_up_ordered)) # num category for plot

      # Guardar dataset
      save(bar_data_up_ordered, file = paste0(outdir, "UP_GO_", plot_name, ".RData"))

      # agregar colores para la grafica
      bar_data_up_ordered_mod <- left_join(bar_data_up_ordered, Category_colors, by= "category")

      ### ---- UP genes (barplot) ----
      # Generar la grafica
      g.up <- ggplot(bar_data_up_ordered_mod, aes(p.val, reorder(term, -num), fill = category)) +
        geom_bar(stat = "identity") +
        geom_text(
          aes(label = p.val),
          color = "black",
          hjust = 0,
          size = 2.2,
          position = position_dodge(0)
        ) +
        labs(x = "-log10(p-value)", y = NULL) +
        scale_fill_manual(
          name = 'Category',
          labels = unique(bar_data_up_ordered_mod$label),
          values = unique(bar_data_up_ordered_mod$colors)
        ) +
        theme(
          legend.position = "right",
          axis.title.y = element_blank(),
          strip.text.x = element_text(size = 11, face = "bold"),
          strip.background = element_blank()
        ) + theme_classic()

      # Guardar la figura
      ggsave(paste0(figdir, "barplotUP_GO_", plot_name, ".png"),
             plot = g.up + theme_classic(), dpi = 600, width = 10, height = 5)
    } else {
      message("No hay términos GO significativos para genes up-regulados en ", plot_name)
    }
  } else {
    message("No se encontraron resultados significativos para ", plot_name)
  }
}
