#!/usr/bin/env bash
#SBATCH --job-name=uniprot_blast
#SBATCH --partition=pibu_el8
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/uniprot_blast_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/uniprot_blast_%J.err

# General paths
WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
ANNODIR="$WORKDIR/results/final"
ORTHODIR="$WORKDIR/results/orthology"

# Create output directory
mkdir -p "$ORTHODIR"

# Input files from previous filtering step
PROTEIN_FASTA="$ANNODIR/assembly_primary.all.maker.proteins.fasta.renamed.filtered.fasta"
GFF_FILE="$ANNODIR/filtered.genes.renamed.gff3"

# UniProt database
UNIPROT_DB="$COURSEDIR/data/uniprot/uniprot_viridiplantae_reviewed.fa"

# MAKER bin
MAKERBIN="$COURSEDIR/softwares/Maker_v3.01.03/src/bin"

# Load required modules
module load BLAST+/2.15.0-gompi-2021a

# Change to orthology directory
cd "$ORTHODIR" || exit 1

# Copy input files
cp "$PROTEIN_FASTA" maker_proteins.fasta
cp "$GFF_FILE" filtered.maker.filtered.gff3

echo "Starting BLASTP against UniProt database..."

# Run BLASTP against UniProt
blastp -query maker_proteins.fasta \
    -db "$UNIPROT_DB" \
    -num_threads 10 \
    -outfmt 6 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -out blastp_uniprot.outfmt6

echo "Sorting BLAST results to get best hits..."

# Sort to keep only the best hit per query sequence
sort -k1,1 -k12,12gr blastp_uniprot.outfmt6 | sort -u -k1,1 --merge > blastp_uniprot.outfmt6.besthits

echo "Creating backup copies..."

# Create backup copies
cp maker_proteins.fasta maker_proteins.fasta.Uniprot
cp filtered.maker.filtered.gff3 filtered.maker.filtered.gff3.Uniprot

echo "Adding functional annotations to FASTA file..."

# Add functional annotations to FASTA
$MAKERBIN/maker_functional_fasta "$UNIPROT_DB" \
    blastp_uniprot.outfmt6.besthits \
    maker_proteins.fasta > maker_proteins.fasta.Uniprot.functional

echo "Adding functional annotations to GFF3 file..."

# Add functional annotations to GFF3
$MAKERBIN/maker_functional_gff "$UNIPROT_DB" \
    blastp_uniprot.outfmt6.besthits \
    filtered.maker.filtered.gff3 > filtered.maker.filtered.gff3.Uniprot.functional

echo "Generating summary statistics..."

# Generate summary statistics
echo "=== UniProt BLAST Summary ===" > uniprot_blast_summary.txt
echo "Total proteins queried: $(grep -c ">" maker_proteins.fasta)" >> uniprot_blast_summary.txt
echo "Proteins with UniProt hits: $(cut -f1 blastp_uniprot.outfmt6.besthits | sort -u | wc -l)" >> uniprot_blast_summary.txt
echo "Total BLAST hits: $(wc -l < blastp_uniprot.outfmt6)" >> uniprot_blast_summary.txt
echo "" >> uniprot_blast_summary.txt

# Check for characterized vs uncharacterized proteins
echo "Analyzing characterized vs uncharacterized proteins..." >> uniprot_blast_summary.txt
grep -i "uncharacterized" blastp_uniprot.outfmt6.besthits | wc -l > temp_unchar
grep -vi "uncharacterized" blastp_uniprot.outfmt6.besthits | wc -l > temp_char
echo "Hits to characterized proteins: $(cat temp_char)" >> uniprot_blast_summary.txt
echo "Hits to uncharacterized proteins: $(cat temp_unchar)" >> uniprot_blast_summary.txt
rm temp_unchar temp_char

cat uniprot_blast_summary.txt

echo "UniProt BLAST analysis complete!"