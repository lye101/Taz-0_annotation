#!/usr/bin/env bash
#SBATCH --job-name=busco_plots
#SBATCH --partition=pibu_el8
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/busco_plots_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/busco_plots_%J.err

# General paths
WORKDIR="/data/users/${USER}/A_thaliana_Taz-0_organization_and_annotation"
OUTDIR="$WORKDIR/results/busco"
LOGDIR="$WORKDIR/log"

# Create directories
mkdir -p "$LOGDIR"
mkdir -p "$OUTDIR"

cd "$OUTDIR"

# Load BUSCO module
module load BUSCO/5.4.2-foss-2021a

# Generate plots for individual runs
generate_plot.py -wd "$OUTDIR/busco_proteins/run_brassicales_odb10"
generate_plot.py -wd "$OUTDIR/busco_transcripts/run_brassicales_odb10"

# Combined plot
mkdir -p "$OUTDIR/combined_summaries"
cp "$OUTDIR/busco_proteins/run_brassicales_odb10"/short_summary*.txt "$OUTDIR/combined_summaries/"
cp "$OUTDIR/busco_transcripts/run_brassicales_odb10"/short_summary*.txt "$OUTDIR/combined_summaries/"
generate_plot.py -wd "$OUTDIR/combined_summaries"
