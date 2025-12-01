#!/bin/bash

#SBATCH --job-name=parseRM
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/parseRM_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/parseRM_%J.err

# Setting the constant for the directories and required files
WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
INPUTDIR="$WORKDIR/results/01a_EDTA_annotation"
GENOME="assembly.fa"
RMOUT="$WORKDIR/results/01a_EDTA_annotation/assembly.fasta.mod.EDTA.anno/assembly.fasta.mod.out"
PARSER="/data/courses/assembly-annotation-course/CDS_annotation/scripts/05-parseRM.pl"
OUTDIR="$WORKDIR/results/parseRM_results"
LOGDIR="$WORKDIR/log"

# Create directories
mkdir -p "$LOGDIR" "$OUTDIR"

# Load the modules
module load BioPerl/1.7.8-GCCcore-10.3.0

# Run the parser (from the input directory where the .mod.out file is)
cd "$INPUTDIR" || exit 1
perl "$PARSER" -i "$RMOUT" -l 50,1 -v

# Move the result files to output directory
# mv ${GENOME}.mod.EDTA.anno/${GENOME}.mod.out.landscape.*.tab "$OUTDIR/" 2>/dev/null
mv assembly.fasta.mod.EDTA.anno/assembly.fasta.mod.out.landscape.*.tab "$OUTDIR/" 2>/dev/null
echo "Results in: $OUTDIR"