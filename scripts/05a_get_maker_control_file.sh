#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=20
#SBATCH --job-name=gene_annotation
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/edta_%j.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/edta_%j.err


WORKDIR=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/results/02_gene_annotations
mkdir -p $WORKDIR
cd $WORKDIR


apptainer exec --bind $WORKDIR \
    /data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif \
    maker -CTL