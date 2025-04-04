#######
# Script: Carga y procesamiento de datos de expresión génica a partir de alineamiento STAR
# Author: Paola Godoy, Yael Hernández y Ariadna Badía
# Date: 10/03/2025
# Description: El siguiente script carga los datos de conteo de genes generados por STAR (archivos ReadsPerGene.out.tab)
# y los combina en una matriz de cuentas. Además, carga los metadatos asociados a las muestras y guarda los resultados
# en formato CSV y RData para su posterior análisis de expresión diferencial.
# Arguments:
#   - Input: metadata.csv, archivos de conteo de STAR (terminación ReadsPerGene.out.tab)
#   - Output: Matriz de cuentas (CSV y RData)
#######

# --- Load data -----
# Cargar archivos
indir <- "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/STAR_output"
outdir <- "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results"

# Verificar que el directorio existe
if (!dir.exists(indir)) {
  stop("El directorio ", indir, " no existe.")
}

# Opcion B - sin movernos de carpeta
files <- dir(indir, pattern = "ReadsPerGene.out.tab")

# Verificar que hay archivos
if (length(files) == 0) {
  stop("No se encontraron archivos 'ReadsPerGene.out.tab' en ", indir)
}

# Crear matriz de cuentas
counts <- c() # esta sera la matriz
for(i in seq_along(files)){
  file_path <- file.path(indir, files[i])
  if (!file.exists(file_path)) {
    stop("El archivo ", file_path, " no existe.")
  }
  x <- read.table(file = file_path, sep = "\t", header = F, as.is = T)
  # as.is para no convertir tipo de datos
  counts <- cbind(counts, x[,2])
}

# Convertir a formato dataframe
counts <- as.data.frame(counts)
rownames(counts) <- x[,1] # Renombrar las filas con el nombre de los genes

# Asignar nombres de muestras a la matriz de cuentas
sample_names <- sub("_ReadsPerGene.out.tab", "", files)
colnames(counts) <- sample_names

# Eliminar las 4 primeras filas
counts <- counts[-c(1:4),]

# Cargar Metadatos
metadata <- read.csv("/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/data/metadata.csv", header = TRUE)

# Almacenar metadata y matriz de cuentas
save(metadata, counts, file = paste0(outdir, "/counts/raw_counts.RData"))
write.csv(counts, file = paste0(outdir, "/counts/raw_counts.csv"))

# Guardar informacion de ejecucion
sessionInfo()
