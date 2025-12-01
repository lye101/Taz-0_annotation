#!/usr/bin/env bash
#SBATCH --job-name=maker_merger
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/maker_merger_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/maker_merger_%J.err

#General  path
WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
LOGDIR="$WORKDIR/log"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"

#path to the directory of the input files
DATASTORE="$WORKDIR/results/02_gene_annotations/assembly.maker.output/assembly_datastore"
MASTERINDEXFILE="$WORKDIR/results/02_gene_annotations/assembly.maker.output/assembly_master_datastore_index.log"

#Path to the program used to merge
MAKERBIN="$COURSEDIR/softwares/Maker_v3.01.03/src/bin"
# Create output directory if it doesn't exist
mkdir -p "$WORKDIR/results/merge_result"

# Merge GFF with sequences
$MAKERBIN/gff3_merge -s -d $MASTERINDEXFILE > "$WORKDIR/results/merge_result/assembly.all.maker.gff"

# Merge GFF without sequences
$MAKERBIN/gff3_merge -n -s -d $MASTERINDEXFILE > "$WORKDIR/results/merge_result/assembly.all.maker.noseq.gff"

# Merge FASTA files
$MAKERBIN/fasta_merge -d $MASTERINDEXFILE -o "$WORKDIR/results/merge_result/assembly_primary"
