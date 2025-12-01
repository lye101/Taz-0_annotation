#!/usr/bin/env bash
#SBATCH --job-name=busco
#SBATCH --partition=pibu_el8
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/busco_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/busco_%J.err

# BUSCO: Quality Assessment of Gene Annotations

WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
ANNODIR="$WORKDIR/results/final"
BUSCODIR="$WORKDIR/results/busco"

# Create BUSCO directory
mkdir -p "$BUSCODIR"

cd "$BUSCODIR" || exit 1

# Load BUSCO module
module load BUSCO/5.4.2-foss-2021a

# Input files
PROTEIN_FASTA="$ANNODIR/maker_proteins.longest.fasta"
TRANSCRIPT_FASTA="$ANNODIR/maker_transcripts.longest.fasta"

echo "Running BUSCO on protein sequences..."

# Run BUSCO on proteins
busco -i "$PROTEIN_FASTA" \
    -l brassicales_odb10 \
    -o busco_proteins \
    -m proteins \
    -c 8

echo "Running BUSCO on transcript sequences..."

# Run BUSCO on transcripts
busco -i "$TRANSCRIPT_FASTA" \
    -l brassicales_odb10 \
    -o busco_transcripts \
    -m transcriptome \
    -c 8

echo ""
echo "BUSCO analysis complete!"
echo "Results are in: $BUSCODIR"
echo ""
echo "Protein BUSCO summary:"
cat busco_proteins/short_summary.specific.brassicales_odb10.busco_proteins.txt
echo ""
echo "Transcript BUSCO summary:"
cat busco_transcripts/short_summary.specific.brassicales_odb10.busco_transcripts.txt