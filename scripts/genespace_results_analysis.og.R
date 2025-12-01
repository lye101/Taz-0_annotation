#!/usr/bin/env Rscript

# GENESPACE Results Analysis
# Analyze core, accessory, and unique orthogroups and genes

library(GENESPACE)
library(data.table)
library(ggplot2)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("Usage: Rscript analyze_genespace_results.R <genespace_working_directory>")
}

wd <- args[1]
pangenome_file <- file.path(wd, "pangenome_matrix.rds")

# Check if pangenome file exists
if (!file.exists(pangenome_file)) {
    stop("Error: pangenome_matrix.rds not found in ", wd)
}

cat("Loading pangenome matrix...\n")
pangenome <- readRDS(pangenome_file)

cat("Analyzing orthogroups...\n\n")

# Get genome names (excluding TAIR10 reference)
genomes <- names(pangenome)
cat("Genomes analyzed:", paste(genomes, collapse = ", "), "\n\n")

# Extract the orthogroup matrix
# Each row is an orthogroup, each column is a genome
# Values are gene counts (0 = absent, >0 = present)

# Calculate presence/absence across genomes
ortho_matrix <- pangenome

# Convert to presence/absence (1 if genes present, 0 if absent)
presence_matrix <- lapply(ortho_matrix, function(x) as.numeric(x > 0))
presence_df <- as.data.frame(presence_matrix)

# Calculate number of genomes each orthogroup is present in
num_genomes <- rowSums(presence_df)
total_genomes <- ncol(presence_df)

cat("Total number of orthogroups:", nrow(presence_df), "\n")
cat("Total number of genomes:", total_genomes, "\n\n")

# Classify orthogroups
core_orthogroups <- num_genomes == total_genomes
accessory_orthogroups <- num_genomes > 1 & num_genomes < total_genomes
unique_orthogroups <- num_genomes == 1

cat("=== Orthogroup Classification ===\n")
cat("Core orthogroups (present in all genomes):", sum(core_orthogroups), "\n")
cat("Accessory orthogroups (present in some genomes):", sum(accessory_orthogroups), "\n")
cat("Unique orthogroups (present in only one genome):", sum(unique_orthogroups), "\n\n")

# Count genes per genome
cat("=== Genes per Genome ===\n")
gene_counts <- sapply(ortho_matrix, function(x) sum(x, na.rm = TRUE))
for (genome in names(gene_counts)) {
    cat(sprintf("  %s: %d genes\n", genome, gene_counts[genome]))
}
cat("\n")

# Core genes per genome
cat("=== Core Genes per Genome ===\n")
core_gene_counts <- sapply(ortho_matrix[core_orthogroups, ], function(x) sum(x, na.rm = TRUE))
for (genome in names(core_gene_counts)) {
    cat(sprintf("  %s: %d core genes\n", genome, core_gene_counts[genome]))
}
cat("\n")

# Accessory genes per genome
cat("=== Accessory Genes per Genome ===\n")
accessory_gene_counts <- sapply(ortho_matrix[accessory_orthogroups, ], function(x) sum(x, na.rm = TRUE))
for (genome in names(accessory_gene_counts)) {
    cat(sprintf("  %s: %d accessory genes\n", genome, accessory_gene_counts[genome]))
}
cat("\n")

# Unique genes per genome
cat("=== Unique Genes per Genome ===\n")
unique_gene_counts <- sapply(ortho_matrix[unique_orthogroups, ], function(x) sum(x, na.rm = TRUE))
for (genome in names(unique_gene_counts)) {
    cat(sprintf("  %s: %d unique genes\n", genome, unique_gene_counts[genome]))
}
cat("\n")

# Create summary table
summary_table <- data.frame(
    Genome = names(gene_counts),
    Total_Genes = as.numeric(gene_counts),
    Core_Genes = as.numeric(core_gene_counts),
    Accessory_Genes = as.numeric(accessory_gene_counts),
    Unique_Genes = as.numeric(unique_gene_counts),
    Core_Percent = round(100 * core_gene_counts / gene_counts, 2),
    Accessory_Percent = round(100 * accessory_gene_counts / gene_counts, 2),
    Unique_Percent = round(100 * unique_gene_counts / gene_counts, 2)
)

cat("=== Summary Table ===\n")
print(summary_table)
cat("\n")

# Save summary table
write.csv(summary_table, file.path(wd, "pangenome_summary.csv"), row.names = FALSE)
cat("Summary table saved to:", file.path(wd, "pangenome_summary.csv"), "\n\n")

# Create visualization
cat("Creating visualization...\n")

# Reshape data for plotting
plot_data <- data.frame(
    Genome = rep(summary_table$Genome, 3),
    Category = rep(c("Core", "Accessory", "Unique"), each = nrow(summary_table)),
    Count = c(summary_table$Core_Genes, summary_table$Accessory_Genes, summary_table$Unique_Genes)
)

plot_data$Category <- factor(plot_data$Category, levels = c("Core", "Accessory", "Unique"))

# Create stacked bar plot
p <- ggplot(plot_data, aes(x = Genome, y = Count, fill = Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Core" = "#2E86AB", "Accessory" = "#A23B72", "Unique" = "#F18F01")) +
    labs(
        title = "Pan-genome Composition",
        x = "Genome",
        y = "Number of Genes",
        fill = "Gene Category"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

ggsave(file.path(wd, "pangenome_composition.pdf"), p, width = 8, height = 6)
cat("Visualization saved to:", file.path(wd, "pangenome_composition.pdf"), "\n\n")

# Create presence/absence heatmap data
cat("Creating presence/absence matrix...\n")
write.csv(presence_df, file.path(wd, "orthogroup_presence_absence.csv"), row.names = TRUE)
cat("Presence/absence matrix saved to:", file.path(wd, "orthogroup_presence_absence.csv"), "\n\n")

cat("Analysis complete!\n")