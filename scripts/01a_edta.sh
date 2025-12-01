#!/bin/bash

#SBATCH --time=7-00:00:00
#SBATCH --mem=400GB
#SBATCH --cpus-per-task=20
#SBATCH --job-name=EDTA_lysander
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/edta_%j.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/edta_%j.err

WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
mkdir -p "$WORKDIR/results/01a_EDTA_annotation"
cd "$WORKDIR/results/01a_EDTA_annotation"

apptainer exec \
    --bind /data \
    /data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif \
    EDTA.pl \
    --genome /data/users/ltucker/A_thaliana_Taz-0_assembly/output/assemblies/LJA/assembly.fasta \
    --species others \
    --step all \
    --sensitive 1 \
    --cds "/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/input_data/TAIR10_cds_20110103_representative_gene_model_updated" \
    --anno 1 \
    --threads 20 \
    --force 1