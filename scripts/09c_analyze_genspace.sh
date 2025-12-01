#!/usr/bin/env bash

#SBATCH --job-name=3genespace_analysis
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/genespace_analysis_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/genespace_analysis_%J.err

# General paths
WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
GENESPACE_DIR="$WORKDIR/results/genespace"

# Path to R script
RSCRIPT="$WORKDIR/scripts/genespace_results_analysis.R"

# Make sure R script exists
if [ ! -f "$RSCRIPT" ]; then
    echo "Error: R script not found at $RSCRIPT"
    echo "Please copy genespace_results_analysis.R to $WORKDIR/scripts/"
    exit 1
fi

# Check if pangenome_matrix.rds exists
if [ ! -f "$GENESPACE_DIR/pangenome_matrix.rds" ]; then
    echo "Error: pangenome_matrix.rds not found in $GENESPACE_DIR"
    echo "Please run GENESPACE first (08b_run_genespace.sh)"
    exit 1
fi

echo "Analyzing GENESPACE results..."
echo "Working directory: $GENESPACE_DIR"

# Run analysis using Apptainer
apptainer exec \
    --bind $COURSEDIR \
    --bind $WORKDIR \
    --bind $SCRATCH:/temp \
    $COURSEDIR/containers/genespace_latest.sif \
    Rscript "$RSCRIPT" "$GENESPACE_DIR" "Taz0"

echo ""
echo "Analysis complete!"
echo ""
echo "Generated files:"
echo "  - pangenome_summary.csv (summary statistics)"
echo "  - pangenome_composition.pdf (visualization)"
echo "  - orthogroup_presence_absence.csv (detailed presence/absence matrix)"
echo ""
echo "Review the results in: $GENESPACE_DIR"