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
#You can edit the script since this line
#
# Your job name
#$ -N trimming
#
# Send an email after the job has finished
#$ -m e
#$ -M sc0bly16@hotmail.com
#
#
# If modules are needed, source modules environment (Do not delete the next line):
. /etc/profile.d/modules.sh
#
# Add any modules you might require:
#
module load trimmomatic/0.33
cd /mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/data

# Bucle para procesar los archivos FASTQ
for i in *_1.fastq.gz; do
    echo "Procesando $i y su pareja ${i%_1.fastq.gz}_2.fastq.gz..."
    trimmomatic PE -threads 8 -phred33 \
        $i \
        "${i%_1.fastq.gz}_2.fastq.gz" \
        /mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/TRIM_results/"${i%_1.fastq.gz}_1_trimmed.fq.gz" \
        /mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/TRIM_results/"${i%_1.fastq.gz}_1_unpaired.fq.gz" \
        /mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/TRIM_results/"${i%_1.fastq.gz}_2_trimmed.fq.gz" \
        /mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/TRIM_results/"${i%_1.fastq.gz}_2_unpaired.fq.gz" \
        ILLUMINACLIP:/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo2/TRIM_results/TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:80
done
