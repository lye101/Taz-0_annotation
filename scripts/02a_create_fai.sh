#!/usr/bin/env bash
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=fai_index
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/fai_index_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/fai_index_%J.err
#SBATCH --partition=pibu_el8


#Setting the constant for the directories and required files 
WORKDIR="/data/users/${USER}/A_thaliana_Taz-0_organization_and_annotation"
OUTDIR="$WORKDIR/results/03a_fai_files"
LOGDIR="$WORKDIR/log"

#Assembly path
HIFIASM="/data/users/ltucker/A_thaliana_Taz-0_assembly/output/assemblies/LJA/assembly.fasta"

#Create the dirs if not present
mkdir -p "$LOGDIR"
mkdir -p "$OUTDIR"

# enter the right dir
cd "$OUTDIR"

#Load samtools module
module load SAMtools/1.13-GCC-10.3.0
samtools faidx "$HIFIASM" #Generate FAI index file
cp "${HIFIASM}.fai" "$OUTDIR/"

awk -F'\t' 'BEGIN{OFS="\t"} {$1="S"$1; print}' "$WORKDIR/results/03a_fai_files/assembly.fasta.fai" > "$WORKDIR/results/03a_fai_files/s.assembly.fasta.fai"