#!/usr/bin/env Rscript

# GENESPACE Analysis Script
# This script runs GENESPACE for comparative genomics and synteny analysis

library(GENESPACE)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if working directory argument is provided
if (length(args) == 0) {
    stop("Usage: Rscript genespace_analysis.R <working_directory>")
}

# Get the folder where the genespace workingDirectory is located
wd <- args[1]

cat("Working directory:", wd, "\n")
cat("Initializing GENESPACE...\n")

# Initialize GENESPACE
# NOTE: path2mcscanx should point to the directory containing MCScanX, not the executable itself
gpar <- init_genespace(
    wd = wd, 
    path2mcscanx = "/data/courses/assembly-annotation-course/CDS_annotation/softwares/MCScanX"
)

cat("Running GENESPACE...\n")

# Run GENESPACE
out <- run_genespace(gpar, overwrite = TRUE)

cat("Querying pan-genome...\n")

# Query pan-genome with TAIR10 as reference
pangenome <- query_pangenes(
    out, 
    bed = NULL, 
    refGenome = "TAIR10", 
    transform = TRUE, 
    showArrayMem = TRUE, 
    showNSOrtho = TRUE, 
    maxMem2Show = Inf
)

cat("Saving pangenome matrix...\n")

# Save pangenome object as rds
saveRDS(pangenome, file = file.path(wd, "pangenome_matrix.rds"))

cat("GENESPACE analysis complete!\n")
cat("Output saved to:", file.path(wd, "pangenome_matrix.rds"), "\n")