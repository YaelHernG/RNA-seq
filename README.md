## Activación Específica de Interferón-γ por el Factor de Transcripción Artificial ATF3 en Células Jurkat: Análisis de RNA-seq y Enriquecimiento Funcional

## Integrantes:
- Paola Albarrán Godoy
- Ariadna Angélica Badía Zamudio
- Yael Daniel Hernández González

**Fecha: 4 de abril del 2025**

**Materia: Bioinformatica**

**4to Semestre**

## Resumen:

Los factores de transcripción artificiales (ATFs) representan una herramienta prometedora para modular la expresión génica con fines terapéuticos. En este reporte, caracterizamos un ATF diseñado para activar la expresión de interferón-gamma (IFN-γ) en células Jurkat humanas mediante análisis de RNA-seq. Utilizando datos de 9 transcriptomas (3 condiciones: ATF3, Aart6, Vehículo), identificamos 135 genes diferencialmente expresados en ATF3 vs Vehículo, incluyendo IFNG, STAT1 y moléculas del complejo HLA. El análisis funcional reveló enriquecimiento en vías inmunes como “cell surface receptor signaling pathway” (p.adj = 1.2e-10) y “defense response to symbiont” (p.adj = 2.8e-5), mientras que el control Aart6 mostró términos no inmunes (ej: “neuronal synaptic plasticity”). Estos resultados demuestran que ATF3 activa específicamente programas transcripcionales asociados a IFN-γ, respaldando su potencial en inmunoterapia contra cáncer e infecciones.

[Reporte final](https://yaelherng.github.io/RNA-seq/)

## Pipeline:
1. Descarga de los datos, mediante ENA browser [download_all_rawData](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/download_all_rawData.sge)
2. Análisis de calidad de los datos crudos [qc1](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/qc1.sge)
3. Trimming mediante la herramienta de Trimmomatic [trimming](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/trimming.sh)
4. Análisis de calidad de los datos trimmeados [qc2](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/qc2.sh) [qc2_job](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/qc2.sge)
5. Indexacion del genoma de referencia con STAR [STAR_index](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/STAR_index.sh)
6. Alineamiento y Conteo de los transcriptomas al genoma de referencia, mediante STAR [align](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/align.sh)
7. Importar de datos a R (archivos de cuentas y metadata) y creación de una matriz de cuentas con todos los transcriptomas [data_R](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/data.R)
8. Creacion de archivo *dds* con *DESeq2*
9. Análisis de Expresion Diferencial Genica (DEG)
10. Normalizacion de los datos
11. Deteccion de Batch effect mediante PCA
12. Obtener los resultados de los contrastes de DEG
    [DEG_analysis](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/DEG_analysis.R)
13. Visualizacion de los datos, mediante Heatmap y Volcanoplot [Visualizacion](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/Visualizacion.R)
14. Analisis de Terminos funcionales (GO terms) [GOterms_analysis](https://github.com/YaelHernG/RNA-seq/blob/main/Scripts/GOterms_analysis.R)
