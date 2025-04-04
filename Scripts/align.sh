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
#$ -N star_alignment_human
#
# Send an email after the job has finished
#$ -m e
#$ -M sc0bly16@hotmail.com
#
# If modules are needed, source modules environment (Do not delete the next line):
. /etc/profile.d/modules.sh
#
# Add any modules you might require:
#
module load star/2.7.11a
#
# Directorios
index=/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/STAR_index
trimmed_dir=/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/TRIM_results
output_dir=/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/STAR_output
#
# Alinear las lecturas recortadas
FILES=$trimmed_dir/*_1_trimmed.fq.gz
for f in $FILES
do
    base=$(basename $f _1_trimmed.fq.gz)
    echo "Procesando muestra: $base"
    STAR --runThreadN 12 \
         --genomeDir $index \
         --readFilesIn $f $trimmed_dir/${base}_2_trimmed.fq.gz \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         --readFilesCommand zcat \
         --outFileNamePrefix $output_dir/${base}_
done
