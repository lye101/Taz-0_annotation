#!/usr/bin/env bash
#SBATCH --job-name=2genespace
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/genespace_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/genespace_%J.err

# General paths
WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
GENESPACE_DIR="$WORKDIR/results/genespace"

# Path to R script
RSCRIPT="$WORKDIR/scripts/genespace_analysis.R"

# Make sure R script exists
if [ ! -f "$RSCRIPT" ]; then
    echo "Error: R script not found at $RSCRIPT"
    echo "Please copy genespace_analysis.R to $WORKDIR/scripts/"
    exit 1
fi

echo "Starting GENESPACE analysis..."
echo "Working directory: $GENESPACE_DIR"
echo "Using container: $COURSEDIR/containers/genespace_latest.sif"

# Run GENESPACE using Apptainer
apptainer exec \
    --bind $COURSEDIR \
    --bind $WORKDIR \
    --bind $SCRATCH:/temp \
    $COURSEDIR/containers/genespace_latest.sif \
    Rscript "$RSCRIPT" "$GENESPACE_DIR"

echo ""
echo "GENESPACE analysis complete!"
echo "Results are in: $GENESPACE_DIR"
echo ""
echo "Generated files:"
ls -lh "$GENESPACE_DIR"/*.rds 2>/dev/null || echo "No .rds files found"
echo ""
echo "Check the results directory for:"
echo "  - pangenome_matrix.rds (main results)"
echo "  - dotplots (.rawHits.pdf, syntenicHits.pdf)"
echo "  - riparian plots (synteny visualization)"
echo ""
echo "Next: Run 09c_analyze_genespace_results.sh to analyze the results"