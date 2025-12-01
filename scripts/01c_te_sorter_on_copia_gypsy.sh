#!/bin/bash

#SBATCH --job-name=extract_TE_families
#SBATCH --partition=pibu_el8
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/ltucker/organization_annotation_course/log/extract_TE_families_%J.out
#SBATCH --error=/data/users/ltucker/organization_annotation_course/log/extract_TE_families_%J.err

# Defining the constant for the path and files
WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
TELIB="${WORKDIR}/results/01a_EDTA_annotation/assembly.fasta.mod.EDTA.TElib.fa"

CONTAINER="/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif"
OUTDIR="$WORKDIR/results/TEsorter_Copia_Gypsy"
LOGDIR="$WORKDIR/log"

# Create directories
mkdir -p "$LOGDIR" "$OUTDIR"

#Load the module for SeqKit 
module load SeqKit/2.6.1

# Extract Copia and Gypsy sequences
seqkit grep -r -p "Copia" "$TELIB" > "$OUTDIR/Copia_sequences.fa"
seqkit grep -r -p "Gypsy" "$TELIB" > "$OUTDIR/Gypsy_sequences.fa"


####################################################################################################################################
# Classify with TEsorter
####################################################################################################################################


# Copia classification
apptainer exec --bind /data "$CONTAINER" \
  TEsorter "$OUTDIR/Copia_sequences.fa" -db rexdb-plant -p "$SLURM_CPUS_PER_TASK" -pre "$OUTDIR/Copia"

# Gypsy classification
apptainer exec --bind /data "$CONTAINER" \
  TEsorter "$OUTDIR/Gypsy_sequences.fa" -db rexdb-plant -p "$SLURM_CPUS_PER_TASK" -pre "$OUTDIR/Gypsy"
