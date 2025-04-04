#!/bin/bash
# Use current working directory
#$ -cwd
#
# Join stdout and stderr
#$ -j y
#
# Run job through bash shell
#$ -S /bin/bash
#
# Your job name
#$ -N import_counts
#
# Send an email after the job has finished
#$ -m e
#$ -M sc0bly16@hotmail.com
#
# If modules are needed, source modules environment (Do not delete the next line):
. /etc/profile.d/modules.sh
#
# Add any modules you might require:
module load r/4.0.2  # Cargar la versión de R necesaria
#
# Directorios
indir="/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/STAR_output"
outdir="/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/results"
#
# Crear directorio de salida si no existe
mkdir -p $outdir/counts
#
# Ejecutar el script de R
Rscript /mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/scripts/data.R
#
# Mensaje de finalización
echo "Job completado con éxito."
