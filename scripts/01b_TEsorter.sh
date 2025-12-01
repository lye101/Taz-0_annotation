#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=20
#SBATCH --job-name=TEsorter
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/edta_%j.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/edta_%j.err


WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
mkdir -p "$WORKDIR/results/01b_TEsorter"
cd "$WORKDIR/results/01b_TEsorter"

# use non-intact fasta file from EDTA output in main outdir
ltr_file="$WORKDIR/results/1_EDTA_annotation/assembly.fasta.mod.EDTA.raw/assembly.fasta.mod.LTR.raw.fa"
container=/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif

apptainer exec --bind /data \
    $container \
    TEsorter \
    $ltr_file \
    -db rexdb-plant