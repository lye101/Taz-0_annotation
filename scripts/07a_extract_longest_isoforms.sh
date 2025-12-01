#!/usr/bin/env bash
#SBATCH --job-name=extract_longest
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/extract_longest_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/extract_longest_%J.err

# Extract longest protein and transcript per gene before running BUSCO
# Gene name is everything before -R, whereas -RA, -RB, -RC etc. are different isoforms

WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
ANNODIR="$WORKDIR/results/final"

cd "$ANNODIR" || exit 1

# Input files
PROTEIN_FASTA="assembly_primary.all.maker.proteins.fasta.renamed.filtered.fasta"
TRANSCRIPT_FASTA="assembly_primary.all.maker.transcripts.fasta.renamed.filtered.fasta"

echo "Extracting longest protein isoform per gene..."

# Extract longest protein per gene
awk 'BEGIN{RS=">"; FS="\n"} 
NR>1 {
    split($1, arr, "-R"); 
    gene = arr[1]; 
    seq = ""; 
    for(i=2; i<=NF; i++) seq = seq $i; 
    len = length(seq); 
    if(len > max_len[gene]) {
        max_len[gene] = len; 
        longest_seq[gene] = ">" $0
    }
} 
END {
    for(gene in longest_seq) print longest_seq[gene]
}' "$PROTEIN_FASTA" > maker_proteins.longest.fasta

echo "Extracting longest transcript isoform per gene..."

# Extract longest transcript per gene
awk 'BEGIN{RS=">"; FS="\n"} 
NR>1 {
    split($1, arr, "-R"); 
    gene = arr[1]; 
    seq = ""; 
    for(i=2; i<=NF; i++) seq = seq $i; 
    len = length(seq); 
    if(len > max_len[gene]) {
        max_len[gene] = len; 
        longest_seq[gene] = ">" $0
    }
} 
END {
    for(gene in longest_seq) print longest_seq[gene]
}' "$TRANSCRIPT_FASTA" > maker_transcripts.longest.fasta

echo "Complete!"
echo "Total proteins: $(grep -c ">" "$PROTEIN_FASTA")"
echo "Longest proteins: $(grep -c ">" maker_proteins.longest.fasta)"
echo "Total transcripts: $(grep -c ">" "$TRANSCRIPT_FASTA")"
echo "Longest transcripts: $(grep -c ">" maker_transcripts.longest.fasta)"