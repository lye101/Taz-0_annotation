#!/usr/bin/env bash
#SBATCH --job-name=tair10_blast
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=15G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/tair10_blast_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/tair10_blast_%J.err

# General paths
WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
ORTHODIR="$WORKDIR/results/orthology"

# TAIR10 database
TAIR10_PEP="$COURSEDIR/data/TAIR10_pep_20110103_representative_gene_model"

# Input file
PROTEIN_FASTA="$ORTHODIR/maker_proteins.fasta"

# Load required modules
module load BLAST+/2.15.0-gompi-2021a

# Change to orthology directory
cd "$ORTHODIR" || exit 1

echo "Starting BLASTP against TAIR10 database..."

# Run BLASTP against TAIR10
blastp -query "$PROTEIN_FASTA" \
    -db "$TAIR10_PEP" \
    -num_threads 10 \
    -outfmt 6 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -out blastp_tair10.outfmt6

echo "Sorting BLAST results to get best hits..."

# Sort to keep only the best hit per query sequence
sort -k1,1 -k12,12gr blastp_tair10.outfmt6 | sort -u -k1,1 --merge > blastp_tair10.outfmt6.besthits

echo "Generating TAIR10 BLAST summary..."

# Generate summary statistics
echo "=== TAIR10 BLAST Summary ===" > tair10_blast_summary.txt
echo "Total proteins queried: $(grep -c ">" "$PROTEIN_FASTA")" >> tair10_blast_summary.txt
echo "Proteins with TAIR10 hits: $(cut -f1 blastp_tair10.outfmt6.besthits | sort -u | wc -l)" >> tair10_blast_summary.txt
echo "Total BLAST hits: $(wc -l < blastp_tair10.outfmt6)" >> tair10_blast_summary.txt
echo "" >> tair10_blast_summary.txt

# Check for FLC gene (AT5G10140)
echo "=== FLC Gene (AT5G10140) Search ===" >> tair10_blast_summary.txt
if grep -q "AT5G10140" blastp_tair10.outfmt6.besthits; then
    echo "FLC gene found!" >> tair10_blast_summary.txt
    grep "AT5G10140" blastp_tair10.outfmt6.besthits >> tair10_blast_summary.txt
else
    echo "FLC gene (AT5G10140) not found in best hits" >> tair10_blast_summary.txt
    # Check all hits
    if grep -q "AT5G10140" blastp_tair10.outfmt6; then
        echo "FLC found in all hits (not best hit):" >> tair10_blast_summary.txt
        grep "AT5G10140" blastp_tair10.outfmt6 | head -5 >> tair10_blast_summary.txt
    fi
fi

cat tair10_blast_summary.txt

echo "TAIR10 BLAST analysis complete!"