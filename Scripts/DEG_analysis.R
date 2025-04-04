######## Script: Análisis de expresión diferencial con DESeq2
# Author: Paola Godoy, Yael Hernández y Ariadna Badía
# Date: 13/03/2025
# Description: Este script realiza un análisis de expresión diferencial (DEG) utilizando DESeq2.
# Carga los datos de conteo de genes y metadatos, normaliza los datos, detecta efectos de batch
# y realiza contrastes entre diferentes condiciones experimentales (control, control_isotopo y tratamiento).
# Arguments:
#   - Input: Matriz de cuentas (raw_counts.RData) y metadatos (metadata.csv).
#   - Output: Resultados de expresión diferencial (archivos CSV), gráficos de PCA y archivos RData con los objetos dds y vsd.
#######

# --- Load packages ----------
library(DESeq2)

# --- Load data -----
outdir <- "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/"
figdir <- '/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/figures/'

# Cargar variable "counts", proveniente del script "load_data_inR.R"
load("/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results/counts/raw_counts.RData")

# Eliminar espacios en blanco en los nombres de las muestras y condiciones
metadata$sample_id <- trimws(metadata$sample_id)
metadata$type <- as.factor(metadata$type) # convertir a factor

samples <- metadata$sample_id # Extraer los nombres de los Transcriptomas
metadata$type <- factor(trimws(as.character(metadata$type)))

# --- DEG ----
counts <- counts[which(rowSums(counts) > 10),] # Seleccionamos genes con más de 10 cuentas

# Convertir al formato dds
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~type) # Se hace un DESeqDataSet para realizar un análisis
dim(dds) # checar las dimensiones

##  -- Asignar la referencia y generar contrastes -----
dds$type <- relevel(dds$type, ref = "vehiculo")

## --- Obtener archivo dds ----
dds <- DESeq(dds)


## --- Normalizacion de los datos ---------
# Opcion 2. regularized logarithm or rlog
dds <- DESeq(dds)

# Obtener la lista de coeficientes o contrastes
resultsNames(dds)

## --- Normalizacion de los datos ---------
# Opcion 2. regularized logarithm or rlog
ddslog <- rlog(dds, blind = F)

# Opcion 3. vsd (más rápido para muchas muestras)
vsdata <- vst(dds, blind = F)

## --- Deteccion de batch effect ----
png(file = paste0(figdir, "PCA_rlog.png"))
plt <- plotPCA(ddslog, intgroup = "type")
print(plt)
dev.off()

png(file = paste0(figdir, "PCA_vsd.png"))
plt <- plotPCA(vsdata, intgroup = "type")
print(plt)
dev.off()

## ---- Obtener informacion del contraste 1
res_Aart6_vs_vehiculo <- results(dds, contrast = c("type", "Aart6", "vehiculo"))
summary(res_Aart6_vs_vehiculo)
write.csv(res_Aart6_vs_vehiculo, file = paste0(outdir, 'DE_Aart6_vs_vehiculo.csv'))

## ---- Obtener informacion del contraste 2 ----
res_ATF3_vs_vehiculo <- results(dds, contrast = c("type", "ATF3", "vehiculo"))
summary(res_ATF3_vs_vehiculo)
write.csv(res_ATF3_vs_vehiculo, file = paste0(outdir,'DE_ATF3_vs_vehiculo.csv'))

## ---- Obtener informacion del contraste 3 ----
res_Aart6_vs_ATF3 <- results(dds, contrast = c("type", "Aart6", "ATF3"))
summary(res_Aart6_vs_ATF3)
write.csv(res_Aart6_vs_ATF3, file = paste0(outdir, 'DE_Aart6_vs_ATF3.csv'))

## ---- Obtener informacion del contraste 4 ----
# Usar el nombre exacto del contraste
res_ATF3_vs_Aart6 <- results(dds, contrast = c("type", "ATF3", "Aart6"))
summary(res_ATF3_vs_Aart6)
write.csv(res_ATF3_vs_Aart6, file = paste0(outdir, 'DE_ATF3_vs_Aart6.csv'))

# Guardar la salida del diseno
save(metadata, dds, file = paste0(outdir, 'dds_Times_vs_vehiculo.RData'))
save(metadata, vsdata, file = paste0(outdir, 'vst_Times_vs_vehiculo.RData'))
