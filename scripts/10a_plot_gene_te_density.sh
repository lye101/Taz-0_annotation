#!/usr/bin/env bash
#SBATCH --job-name=visualize_gene_te
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/visualize_gene_te_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/visualize_gene_te_%J.err

# Visualize Gene and TE Annotations Together
# Answers guiding questions from the course

WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
VIZDIR="$WORKDIR/results/visualization"

# Create visualization directory
mkdir -p "$VIZDIR"

# Input files
GENE_GFF="$WORKDIR/results/final/filtered.genes.renamed.gff3"
TE_GFF="$WORKDIR/results/01a_EDTA_annotation/assembly.fasta.mod.EDTA.TEanno.gff3"
GENOME_FAI="$WORKDIR/results/03a_fai_files/assembly.fasta.fai"

# R script
RSCRIPT="$WORKDIR/scripts/visualize_gene_te_density.R"

# Check if input files exist
if [ ! -f "$GENE_GFF" ]; then
    echo "ERROR: Gene GFF file not found: $GENE_GFF"
    exit 1
fi

if [ ! -f "$TE_GFF" ]; then
    echo "ERROR: TE GFF file not found: $TE_GFF"
    exit 1
fi

if [ ! -f "$GENOME_FAI" ]; then
    echo "ERROR: Genome FAI file not found: $GENOME_FAI"
    echo "Creating FAI file..."
    module load SAMtools/1.13-GCC-10.3.0
    samtools faidx "$WORKDIR/results/01a_EDTA_annotation/assembly.fasta.mod"
    GENOME_FAI="$WORKDIR/results/01a_EDTA_annotation/assembly.fasta.mod.fai"
fi

echo "=== Gene and TE Density Visualization ==="
echo "Gene GFF: $GENE_GFF"
echo "TE GFF: $TE_GFF"
echo "Genome FAI: $GENOME_FAI"
echo ""

# Change to visualization directory
cd "$VIZDIR" || exit 1

# Run visualization using Apptainer with R
apptainer exec \
    --bind $WORKDIR \
    --bind $COURSEDIR \
    $COURSEDIR/containers/genespace_latest.sif \
    Rscript "$RSCRIPT" "$GENE_GFF" "$TE_GFF" "$GENOME_FAI" 100000

echo ""
echo "=== Visualization Complete ==="
echo "Output files in: $VIZDIR"
ls -lh "$VIZDIR"/*.pdf
ls -lh "$VIZDIR"/*.csv

echo ""
echo "Generated visualizations:"
echo "  1. gene_te_density_chromosomes.pdf - Density along each chromosome"
echo "  2. gene_vs_te_density_scatter.pdf - Gene vs TE density correlation"
echo "  3. gene_te_average_density.pdf - Overall comparison"
echo "  4. gene_te_density_heatmap.pdf - Heatmap view"
echo ""
echo "Data files:"
echo "  - gene_te_density_summary.csv - Summary statistics"
echo "  - gene_te_density_by_window.csv - Detailed window data"