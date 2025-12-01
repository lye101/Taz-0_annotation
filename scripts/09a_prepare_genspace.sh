#!/usr/bin/env bash
#SBATCH --job-name=genespace_prep
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/genespace_prep_%J.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/genespace_prep_%J.err

set -e

# General paths
WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
ANNODIR="$WORKDIR/results/final"
GENESPACE_DIR="$WORKDIR/results/genespace"

# Your accession name
ACCESSION="Taz0"

# Input files
GFF_FILE="$ANNODIR/filtered.genes.renamed.gff3"
PROTEIN_FASTA="$ANNODIR/assembly_primary.all.maker.proteins.fasta.renamed.filtered.fasta"

# TAIR10 reference files
TAIR10_BED="$COURSEDIR/data/TAIR10.bed"
TAIR10_FA="$COURSEDIR/data/TAIR10.fa"

# Lian et al. data paths
LIAN_GFF="$COURSEDIR/data/Lian_et_al/gene_gff/selected"
LIAN_PROTEIN="$COURSEDIR/data/Lian_et_al/protein/selected"

# Other accessions (same as original script)
OTHER_ACCESSIONS=("Altai_5" "Are_6" "Etna_2")

# Create output directories
mkdir -p "${GENESPACE_DIR}/peptide"
mkdir -p "${GENESPACE_DIR}/bed"

# Process your accession (Taz0)
echo "Processing ${ACCESSION}..."
grep -P "\tgene\t" "${GFF_FILE}" | \
awk 'BEGIN{OFS="\t"} {
    split($9, a, ";");
    split(a[1], b, "=");
    gene_id = b[2];
    gsub(/[:.-]/, "_", gene_id);
    print $1, $4-1, $5, gene_id
}' | sort -k1,1 -k2,2n > "${GENESPACE_DIR}/bed/${ACCESSION}.bed"

awk '
/^>/ {
    if(seq != "") {
        print seq;
    }
    id = substr($1, 2);
    sub(/-R.*/, "", id);
    sub(/\.[0-9]+$/, "", id);
    gsub(/[:.-]/, "_", id);
    print ">" id;
    seq = "";
    next;
}
{
    gsub(/[^A-Za-z*]/, "", $0);
    seq = seq $0;
}
END {
    if(seq != "") {
        print seq;
    }
}
' "${PROTEIN_FASTA}" > "${GENESPACE_DIR}/peptide/${ACCESSION}.fa"

# Process TAIR10
echo "Processing TAIR10..."
awk 'BEGIN{OFS="\t"} {
    gene_id = $4;
    sub(/\.[0-9]+$/, "", gene_id);
    gsub(/[:.-]/, "_", gene_id);
    print $1, $2, $3, gene_id;
}' "${TAIR10_BED}" | sort -k1,1 -k2,2n > "${GENESPACE_DIR}/bed/TAIR10.bed"

awk '
/^>/ {
    if(seq != "") {
        print seq;
    }
    id = substr($1, 2);
    sub(/\.[0-9]+$/, "", id);
    gsub(/[:.-]/, "_", id);
    print ">" id;
    seq = "";
    next;
}
{
    gsub(/[^A-Za-z*]/, "", $0);
    seq = seq $0;
}
END {
    if(seq != "") {
        print seq;
    }
}
' "${TAIR10_FA}" > "${GENESPACE_DIR}/peptide/TAIR10.fa"

# Process other accessions from Lian et al.
for ACC in "${OTHER_ACCESSIONS[@]}"; do
    echo "Processing ${ACC}..."
    ACC_DASH="${ACC//_/-}"
    
    GFF="${LIAN_GFF}/${ACC_DASH}.*.gff"
    GFF=$(ls ${GFF} 2>/dev/null | head -n 1)
    
    PROT="${LIAN_PROTEIN}/${ACC_DASH}.protein.faa"
    PROT=$(ls ${PROT} 2>/dev/null | head -n 1)
    
    if [[ -z "${GFF}" ]] || [[ -z "${PROT}" ]]; then
        echo "Warning: Files not found for ${ACC}, skipping..."
        continue
    fi
    
    grep -P "\tgene\t" "${GFF}" | \
    awk 'BEGIN{OFS="\t"} {
        split($9, a, ";");
        split(a[1], b, "=");
        gene_id = b[2];
        gsub(/[:.-]/, "_", gene_id);
        print $1, $4-1, $5, gene_id
    }' | sort -k1,1 -k2,2n > "${GENESPACE_DIR}/bed/${ACC}.bed"
    
    awk '
    /^>/ {
        if(seq != "") {
            print seq;
        }
        id = substr($1, 2);
        sub(/-R.*/, "", id);
        sub(/\.[0-9]+$/, "", id);
        gsub(/[:.-]/, "_", id);
        print ">" id;
        seq = "";
        next;
    }
    {
        gsub(/[^A-Za-z*]/, "", $0);
        seq = seq $0;
    }
    END {
        if(seq != "") {
            print seq;
        }
    }
    ' "${PROT}" > "${GENESPACE_DIR}/peptide/${ACC}.fa"
done

echo "Done! GENESPACE preparation complete."