#!/usr/bin/env bash
#SBATCH --job-name=agat_stats
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/agat_stats_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/agat_stats_%J.err

# AGAT: Annotation Statistics

WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
ANNODIR="$WORKDIR/results/final"
STATSDIR="$WORKDIR/results/annotation_stats"

# Create stats directory
mkdir -p "$STATSDIR"

cd "$STATSDIR" || exit 1

# Input GFF file
GFF_FILE="$ANNODIR/filtered.genes.renamed.gff3"

echo "Running AGAT statistics on gene annotations..."

# Run AGAT using apptainer
apptainer exec \
    --bind $WORKDIR \
    --bind $COURSEDIR \
    $COURSEDIR/containers/agat_1.5.1--pl5321hdfd78af_0.sif \
    agat_sp_statistics.pl -i "$GFF_FILE" -o annotation.stat

echo ""
echo "AGAT statistics complete!"
echo "Results saved to: $STATSDIR/annotation.stat"
echo ""
echo "Summary statistics:"
cat annotation.stat