#!/usr/bin/env bash

#SBATCH --job-name=filter_gene_annotation
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=15G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/filter_gene_annotation_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/filter_gene_annotation_%J.err

# General path
WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
LOGDIR="$WORKDIR/log"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"

# Path to the directory of the input files
ANNODIR="$WORKDIR/results"
INPUTDIR="$ANNODIR/merge_result/"

# Path to the program used to merge
MAKERBIN="$COURSEDIR/softwares/Maker_v3.01.03/src/bin"

# Define prefix for gene IDs
prefix="Taz0"  # CHANGE THIS TO YOUR PREFERRED PREFIX

# Create log directory if it doesn't exist
mkdir -p "$LOGDIR"

# Create final directory
mkdir -p "$ANNODIR/final"

# Define file names (based on your merge_result directory)
protein="assembly_primary.all.maker.proteins.fasta"
transcript="assembly_primary.all.maker.transcripts.fasta"
gff="assembly.all.maker.noseq.gff"

# Copy files from merge_result to final directory
cp "$INPUTDIR/$gff" "$ANNODIR/final/${gff}.renamed.gff"
cp "$INPUTDIR/$protein" "$ANNODIR/final/${protein}.renamed.fasta"
cp "$INPUTDIR/$transcript" "$ANNODIR/final/${transcript}.renamed.fasta"

# Change to final directory
cd "$ANNODIR/final" || exit 1

# 1. To assign clean, consistent IDs to the gene models, use MAKER's ID mapping tools
$MAKERBIN/maker_map_ids --prefix $prefix --justify 7 ${gff}.renamed.gff > id.map
$MAKERBIN/map_gff_ids id.map ${gff}.renamed.gff
$MAKERBIN/map_fasta_ids id.map ${protein}.renamed.fasta
$MAKERBIN/map_fasta_ids id.map ${transcript}.renamed.fasta

# 2. Run InterProScan on the Protein File
apptainer exec \
  --bind $COURSEDIR/data/interproscan-5.70-102.0/data:/opt/interproscan/data \
  --bind $WORKDIR \
  --bind $COURSEDIR \
  --bind $SCRATCH:/temp \
  $COURSEDIR/containers/interproscan_latest.sif \
  /opt/interproscan/interproscan.sh \
  -appl pfam --disable-precalc -f TSV --goterms --iprlookup --seqtype p \
  -i ${protein}.renamed.fasta -o output.iprscan

# 3. Update GFF with InterProScan Results
$MAKERBIN/ipr_update_gff ${gff}.renamed.gff output.iprscan > ${gff}.renamed.iprscan.gff

# 4. Calculate AED Values
perl $MAKERBIN/AED_cdf_generator.pl -b 0.025 ${gff}.renamed.gff > assembly.all.maker.renamed.gff.AED.txt

# 5. Filter the GFF File for Quality
perl $MAKERBIN/quality_filter.pl -s ${gff}.renamed.iprscan.gff > ${gff}_iprscan_quality_filtered.gff

# 6. Filter the GFF File for Gene Features
# We only want to keep gene features in the third column of the gff file
grep -P "\tgene\t|\tCDS\t|\texon\t|\tfive_prime_UTR\t|\tthree_prime_UTR\t|\tmRNA\t" ${gff}_iprscan_quality_filtered.gff > filtered.genes.renamed.gff3

# Check
cut -f3 filtered.genes.renamed.gff3 | sort | uniq

# 7. Extract mRNA Sequences and Filter FASTA Files
module load UCSC-Utils/448-foss-2021a
module load MariaDB/10.6.4-GCC-10.3.0
grep -P "\tmRNA\t" filtered.genes.renamed.gff3 | awk '{print $9}' | cut -d ';' -f1 | sed 's/ID=//g' > list.txt
faSomeRecords ${transcript}.renamed.fasta list.txt ${transcript}.renamed.filtered.fasta
faSomeRecords ${protein}.renamed.fasta list.txt ${protein}.renamed.filtered.fasta