#!/usr/bin/env Rscript

# GENESPACE Pangenome Analysis - Pipeline Version
# Analyzes core, accessory, and species-specific orthogroups and genes
# with focal genome comparisons and synteny filtering

library(data.table)
library(ggplot2)

# -----------------
# Command line arguments
# -----------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    cat("Usage: Rscript genespace_analysis_pipeline.R <genespace_working_directory> [focal_genome]\n")
    cat("\n")
    cat("Arguments:\n")
    cat("  genespace_working_directory : Path to directory containing pangenome_matrix.rds\n")
    cat("  focal_genome (optional)     : Name of focal genome for comparative analysis (default: auto-detect)\n")
    cat("\n")
    cat("Example:\n")
    cat("  Rscript genespace_analysis_pipeline.R ./genespace_results Taz0\n")
    cat("\n")
    stop("No arguments provided")
}

wd <- args[1]
pangenome_file <- file.path(wd, "pangenome_matrix.rds")
focal_genome <- "Taz0"

# Check if pangenome file exists
if (!file.exists(pangenome_file)) {
    stop("Error: pangenome_matrix.rds not found in ", wd)
}

# Check if working directory exists
if (!dir.exists(wd)) {
    stop("Error: Working directory does not exist: ", wd)
}

cat("======================================\n")
cat("GENESPACE Pangenome Analysis\n")
cat("======================================\n\n")
cat("Working directory:", wd, "\n")

# -----------------
# 0) Load data
# -----------------
cat("\nLoading pangenome matrix...\n")
pangenome <- readRDS(pangenome_file)

# genome columns are list-columns produced by query_pangenes()
genome_cols <- names(pangenome)[sapply(pangenome, is.list)]

if (length(genome_cols) == 0) {
    stop("Error: No genome columns (list-columns) found in pangenome_matrix.rds")
}

# Determine focal genome
if (length(args) >= 2) {
    focal_genome <- args[2]
    if (!focal_genome %in% genome_cols) {
        cat("\nWarning: Specified focal genome '", focal_genome, "' not found in genome columns.\n", sep = "")
        cat("Available genomes:", paste(genome_cols, collapse = ", "), "\n")
        stop("Please specify a valid focal genome name")
    }
} else {
    # Auto-detect: use first non-TAIR10 genome or first genome
    if ("TAIR10" %in% genome_cols && length(genome_cols) > 1) {
        focal_genome <- genome_cols[genome_cols != "TAIR10"][1]
        cat("Auto-detected focal genome:", focal_genome, "\n")
    } else {
        focal_genome <- genome_cols[1]
        cat("Using first genome as focal:", focal_genome, "\n")
    }
}

cat("Focal genome:", focal_genome, "\n")
cat("All genomes:", paste(genome_cols, collapse = ", "), "\n\n")

# Convert to data.table
pg <- as.data.table(pangenome)

# -----------------
# 1) Clean data - Remove out-of-synteny genes
# -----------------
cat("Filtering out-of-synteny genes (IDs ending with '*')...\n")

clean_gene_list <- function(v) {
  if (is.null(v) || length(v) == 0) {
    return(character(0))
  }
  v <- as.character(v)
  v <- v[!is.na(v)]
  v <- trimws(v)
  v <- v[!grepl("\\*$", v)] # drop genes ending with '*'
  unique(v)
}

# Apply cleaning to all genome columns
for (col in genome_cols) {
  pg[[col]] <- lapply(pg[[col]], clean_gene_list)
}

# Remove orthogroups that are now empty in ALL genomes after cleaning
pg[, total_genes := Reduce(`+`, lapply(.SD, function(x) sapply(x, length))), 
   .SDcols = genome_cols]
n_orthogroups_before <- nrow(pg)
pg <- pg[total_genes > 0]
pg[, total_genes := NULL]
n_orthogroups_removed <- n_orthogroups_before - nrow(pg)

cat("  - Orthogroups before filtering:", n_orthogroups_before, "\n")
cat("  - Orthogroups removed (empty after filtering):", n_orthogroups_removed, "\n")
cat("  - Orthogroups retained:", nrow(pg), "\n\n")

# -----------------
# 2) Presence/absence per genome (orthogroup-level)
# -----------------
cat("Calculating presence/absence matrix...\n")

presence_tbl <- copy(pg)

# Create presence/absence columns
for (col in genome_cols) {
  presence_tbl[[col]] <- sapply(pg[[col]], function(x) length(x) > 0)
}

# Keep only pgID and presence columns
presence_cols <- c("pgID", genome_cols)
presence_tbl <- presence_tbl[, ..presence_cols]

# -----------------
# 3) Global orthogroup flags
# -----------------
n_genomes <- length(genome_cols)

cat("Classifying orthogroups...\n")

# Calculate number of genomes each orthogroup is present in
presence_tbl[, n_present := Reduce(`+`, lapply(.SD, as.integer)), 
             .SDcols = genome_cols]

# Add flags
presence_tbl[, `:=`(
  # Core = present in ALL genomes
  is_core_all = (n_present == n_genomes),
  
  # Accessory = present in SOME genomes (not all)
  is_accessory = (n_present < n_genomes),
  
  # Is this orthogroup in our focal genome?
  present_in_focal = get(focal_genome),
  
  # Focal-specific = ONLY in focal genome
  focal_specific = (get(focal_genome) & n_present == 1),
  
  # Lost in focal = present in other genomes but NOT focal
  lost_in_focal = (!get(focal_genome) & n_present > 0),
  
  # Almost-core lost in focal = in all genomes EXCEPT focal
  lost_core_without_focal = (!get(focal_genome) & n_present == (n_genomes - 1))
)]

# Check if TAIR10 column exists before creating lost_vs_TAIR10
if ("TAIR10" %in% genome_cols) {
  presence_tbl[, lost_vs_TAIR10 := (!present_in_focal) & TAIR10]
} else {
  presence_tbl[, lost_vs_TAIR10 := FALSE]
}

# Add category labels
presence_tbl[, category := "accessory"]
presence_tbl[is_core_all == TRUE, category := "core"]
presence_tbl[n_present == 1, category := "species_specific"]

pg_flags <- presence_tbl

# Print orthogroup classification summary
cat("\n=== Orthogroup Classification ===\n")
cat("Total orthogroups:", nrow(pg_flags), "\n")
cat("Core orthogroups (all genomes):", sum(pg_flags$is_core_all), "\n")
cat("Accessory orthogroups (some genomes):", sum(pg_flags$is_accessory & pg_flags$n_present > 1), "\n")
cat("Species-specific (1 genome):", sum(pg_flags$n_present == 1), "\n")
cat("\n")

cat("=== Focal Genome ('", focal_genome, "') Analysis ===\n", sep = "")
cat("Orthogroups present in focal:", sum(pg_flags$present_in_focal), "\n")
cat("Focal-specific orthogroups:", sum(pg_flags$focal_specific), "\n")
cat("Orthogroups lost in focal:", sum(pg_flags$lost_in_focal), "\n")
cat("Almost-core lost in focal:", sum(pg_flags$lost_core_without_focal), "\n")
if ("TAIR10" %in% genome_cols) {
  cat("Lost compared to TAIR10:", sum(pg_flags$lost_vs_TAIR10), "\n")
}
cat("\n")

# -----------------
# 4) GENE counts per orthogroup × genome
# -----------------
cat("Counting genes per orthogroup...\n")

count_genes <- function(gene_list) {
  if (is.null(gene_list) || length(gene_list) == 0) {
    return(0L)
  }
  length(unique(gene_list))
}

# Create gene counts matrix
gene_counts_matrix <- copy(pg[, .(pgID)])

# Apply the function to all genome columns
for (col in genome_cols) {
  gene_counts_matrix[[col]] <- sapply(pg[[col]], count_genes)
}

# -----------------
# 5) GENE counts per genome & category
# -----------------
cat("Summarizing gene counts by category...\n")

# Add category to the matrix
gene_counts_w_cat <- merge(
  pg_flags[, .(pgID, category)],
  gene_counts_matrix,
  by = "pgID"
)

# Melt to long format for summarization
gene_counts_long <- melt(
  gene_counts_w_cat,
  id.vars = c("pgID", "category"),
  measure.vars = genome_cols,
  variable.name = "genome",
  value.name = "gene_count"
)

# Calculate total genes per genome by category
gene_by_cat <- gene_counts_long[, .(gene_count = sum(gene_count)), 
                                 by = .(genome, category)]

# Total genes per genome (sum across all categories)
gene_totals <- gene_counts_long[, .(gene_total = sum(gene_count)), 
                                by = genome]

# Spread categories for a per-genome wide table
gene_counts_per_genome <- dcast(
  gene_by_cat,
  genome ~ category,
  value.var = "gene_count",
  fill = 0
)

# Rename columns
setnames(gene_counts_per_genome, 
         old = c("core", "accessory", "species_specific"),
         new = c("gene_core", "gene_accessory", "gene_specific"),
         skip_absent = TRUE)

# Add totals
gene_counts_per_genome <- merge(gene_counts_per_genome, gene_totals, by = "genome")

# Calculate percentages
gene_counts_per_genome[, `:=`(
  percent_core = round(100 * gene_core / pmax(gene_total, 1), 2),
  percent_accessory = round(100 * gene_accessory / pmax(gene_total, 1), 2),
  percent_specific = round(100 * gene_specific / pmax(gene_total, 1), 2)
)]

cat("\n=== Gene Counts per Genome ===\n")
print(gene_counts_per_genome)
cat("\n")

# Save summary table
summary_file <- file.path(wd, "pangenome_gene_summary.csv")
fwrite(gene_counts_per_genome, summary_file)
cat("Gene summary saved to:", summary_file, "\n\n")

# -----------------
# 6) Pangenome frequency plot: OGs and genes in n genomes
# -----------------
cat("Generating frequency distribution data...\n")

# Count orthogroups by number of genomes they're present in
og_freq <- pg_flags[, .(count = .N), by = n_present]
og_freq[, type := "Orthogroups"]

# Count genes: how many present in orthogroups present in n genomes
# Create a long table of (pgID, genome, gene) and attach the orthogroup-level n_present

# First, create long format with all genes
all_genes_list <- list()

for (genome_name in genome_cols) {
  # Extract genes for this genome
  genes_data <- pg[, .(pgID, genes_list = get(genome_name))]
  
  # Expand the list column
  genes_expanded <- genes_data[, .(gene = unlist(genes_list)), by = pgID]
  
  # Filter out NAs and empty strings
  genes_expanded <- genes_expanded[!is.na(gene) & nchar(as.character(gene)) > 0]
  
  # Add genome name
  genes_expanded[, genome := genome_name]
  
  all_genes_list[[genome_name]] <- genes_expanded
}

# Combine all genomes
all_genes_with_genome <- rbindlist(all_genes_list)

# Get unique gene occurrences per pgID and genome
all_genes_unique <- unique(all_genes_with_genome[, .(pgID, genome, gene)])

# Add n_present from pg_flags
all_genes_with_presence <- merge(
  all_genes_unique,
  pg_flags[, .(pgID, n_present)],
  by = "pgID"
)

# Count genes by the n_present of their orthogroup
gene_freq <- unique(all_genes_with_presence[, .(pgID, gene, n_present)])[
  , .(count = .N), by = n_present
]
gene_freq[, type := "Genes"]

# Combine both
freq_data <- rbind(og_freq, gene_freq)

# Create summary table
freq_summary <- dcast(
  freq_data,
  n_present ~ type,
  value.var = "count",
  fill = 0
)

# Add labels
freq_summary[, label := ifelse(
  n_present == n_genomes, 
  "Core (all genomes)",
  ifelse(
    n_present == 1,
    "Species-specific (1 genome)",
    paste0("Shared (", n_present, " genomes)")
  )
)]

# Reorder columns
setcolorder(freq_summary, c("n_present", "label", "Orthogroups", "Genes"))

# Sort by n_present descending
setorder(freq_summary, -n_present)

cat("\n=== Frequency Distribution ===\n")
print(freq_summary)
cat("\n")

# Save frequency summary
freq_file <- file.path(wd, "pangenome_frequency_summary.csv")
fwrite(freq_summary, freq_file)
cat("Frequency summary saved to:", freq_file, "\n\n")

# -----------------
# 7) Create visualizations
# -----------------
cat("Creating visualizations...\n")

# Plot 1: Frequency distribution
p1 <- ggplot(freq_data, aes(x = n_present, y = count, fill = type)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  scale_x_continuous(
    breaks = 1:n_genomes,
    labels = 1:n_genomes
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("Genes" = "#0072B2", "Orthogroups" = "#D55E00")) +
  labs(
    x = "Number of genomes",
    y = "Count",
    fill = NULL,
    title = "Pangenome composition: distribution across genomes",
    subtitle = paste("Total genomes:", n_genomes, "| Focal genome:", focal_genome)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

plot1_file <- file.path(wd, "pangenome_frequency_plot.pdf")
ggsave(plot1_file, plot = p1, width = 10, height = 6)
cat("  - Frequency plot saved to:", plot1_file, "\n")

# Plot 2: Gene composition per genome (stacked bar)
plot_data <- melt(
  gene_counts_per_genome[, .(genome, gene_core, gene_accessory, gene_specific)],
  id.vars = "genome",
  variable.name = "category",
  value.name = "count"
)

# Clean up category names
plot_data[, category := gsub("gene_", "", category)]
plot_data[, category := factor(category, levels = c("core", "accessory", "specific"))]

p2 <- ggplot(plot_data, aes(x = genome, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c("core" = "#2E86AB", "accessory" = "#A23B72", "specific" = "#F18F01"),
    labels = c("Core", "Accessory", "Species-specific")
  ) +
  labs(
    title = "Gene composition per genome",
    subtitle = paste("Focal genome:", focal_genome),
    x = "Genome",
    y = "Number of genes",
    fill = "Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

plot2_file <- file.path(wd, "pangenome_composition_per_genome.pdf")
ggsave(plot2_file, plot = p2, width = 10, height = 6)
cat("  - Composition plot saved to:", plot2_file, "\n")

# -----------------
# 8) Save presence/absence matrix
# -----------------
cat("\nSaving presence/absence matrix...\n")
presence_file <- file.path(wd, "orthogroup_presence_absence.csv")
fwrite(presence_tbl, presence_file)
cat("Presence/absence matrix saved to:", presence_file, "\n")

# -----------------
# Summary
# -----------------
cat("\n======================================\n")
cat("Analysis Complete!\n")
cat("======================================\n\n")
cat("Generated files in:", wd, "\n")
cat("  - pangenome_gene_summary.csv\n")
cat("  - pangenome_frequency_summary.csv\n")
cat("  - pangenome_frequency_plot.pdf\n")
cat("  - pangenome_composition_per_genome.pdf\n")
cat("  - orthogroup_presence_absence.csv\n")
cat("\n")