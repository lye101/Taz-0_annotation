# Organization and Annotation of Eukaryote Genomes 🧬🌱

Resources, code, and materials for the genome annotation of *Arabidopsis thaliana* accession **Taz-0**.

---

## Overview

This repository contains scripts for the annotation of an *Arabidopsis thaliana* Taz-0 assembly. The assembly was generated during the Genome and Transcriptome Assembly course as part of the MSc Bioinformatics and Computational Biology program.

The project aims to perform comprehensive genome annotation of the Taz-0 accession, with particular focus on:

- **Transposable element annotation** to infer evolutionary patterns
- **Pan-genome analysis** comparing Taz-0 with accessions from diverse geographic regions to explore genomic variation across populations

---

## Sample 🌱

**Accession:** Taz-0 (*Arabidopsis thaliana*)

### References

- Lian Q *et al.* (2024) A pan-genome of 69 *Arabidopsis thaliana* accessions reveals a conserved genome structure throughout the global species range. *Nature Genetics* 56:982–991.
- Jiao WB & Schneeberger K (2020) Chromosome-level assemblies of multiple *Arabidopsis* genomes reveal hotspots of rearrangements with altered evolutionary dynamics. *Nature Communications* 11:989.

---

## Workflow 👩🏽‍💻

### 1. Transposable Element Annotation
Annotate TEs using **EDTA** and classify full-length LTR retrotransposons.

### 2. TE Visualization
Generate density and distribution plots across the genome using **circlize** in R.

### 3. TE Classification Refinement
Extract Copia and Gypsy superfamilies; refine clade classification using **TEsorter**.

### 4. TE Dynamics
Estimate insertion ages and divergence from **RepeatMasker** outputs; visualize TE landscapes.

### 5. Gene Annotation
Perform evidence-based gene prediction with **MAKER**, integrating transcriptomic data, protein homology, and *ab initio* predictions. Filter for high-quality gene models.

### 6. Quality Assessment
Evaluate annotation completeness using **BUSCO**; generate summary statistics with **AGAT**.

### 7. Functional Annotation
Assign putative functions via **BLASTP** against UniProt and TAIR10; map annotations to GFF3 and FASTA files.

### 8. Pan-genome Analysis
Compare Taz-0 to other accessions using **GENESPACE** to identify core, accessory, and unique genes; visualize orthogroup distributions.

---

## Tools 🛠️

| Tool | Version |
|------|---------|
| EDTA | 2.2 |
| TEsorter | 1.3.0 |
| SAMtools | 1.13 |
| SeqKit | 2.6.1 |
| MAKER | 3.01.03 |
| OpenMPI | 4.1.1 |
| BioPerl | 1.7.8 |
| AUGUSTUS | 3.4.0 |
| R | 4.5.0 |
| InterProScan | 5.70-102.0 |
| AGAT | 1.5.1 |
| BLAST+ | 2.15.0 |

---

## Repository Structure

```
/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation
├── input_data/
├── log/
├── results/
└── scripts/
```

---
