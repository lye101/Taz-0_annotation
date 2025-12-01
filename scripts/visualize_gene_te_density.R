#!/usr/bin/env Rscript

# Visualize Gene and TE Annotations Together
# Creates a circos plot with multiple tracks showing:
# - TE density
# - Gene density
# - Gene annotations
# - TE annotations

library(circlize)
library(GenomicRanges)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    stop("Usage: Rscript visualize_genes_TEs.R <genome.fasta.fai> <genes.gff3> <TEs.gff3> <output.pdf>")
}

fai_file <- args[1]
genes_gff <- args[2]
tes_gff <- args[3]
output_pdf <- args[4]

cat("Loading data...\n")
cat("Genome index:", fai_file, "\n")
cat("Gene annotations:", genes_gff, "\n")
cat("TE annotations:", tes_gff, "\n")
cat("Output:", output_pdf, "\n\n")

# Read chromosome lengths from FAI file
fai <- read.table(fai_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fai) <- c("chr", "length", "offset", "linebases", "linewidth")

# Filter for main chromosomes (adjust if needed)
# For Arabidopsis: Chr1, Chr2, Chr3, Chr4, Chr5
main_chrs <- fai[grep("^Chr[0-9]", fai$chr), ]
if (nrow(main_chrs) == 0) {
    cat("No 'Chr' chromosomes found, using all sequences > 1Mb\n")
    main_chrs <- fai[fai$length > 1000000, ]
}

cat("Chromosomes to plot:\n")
print(main_chrs[, c("chr", "length")])

# Read gene annotations
cat("\nReading gene annotations...\n")
genes <- read.table(genes_gff, header = FALSE, sep = "\t", 
                    comment.char = "#", stringsAsFactors = FALSE,
                    fill = TRUE)
colnames(genes) <- c("chr", "source", "type", "start", "end", 
                     "score", "strand", "phase", "attributes")

# Filter for gene features
genes <- genes[genes$type == "gene", ]
genes <- genes[genes$chr %in% main_chrs$chr, ]
cat("Number of genes:", nrow(genes), "\n")

# Read TE annotations
cat("Reading TE annotations...\n")
tes <- read.table(tes_gff, header = FALSE, sep = "\t", 
                  comment.char = "#", stringsAsFactors = FALSE,
                  fill = TRUE)
colnames(tes) <- c("chr", "source", "type", "start", "end", 
                   "score", "strand", "phase", "attributes")

# Filter TEs
tes <- tes[tes$chr %in% main_chrs$chr, ]
cat("Number of TEs:", nrow(tes), "\n")

# Calculate densities in windows
cat("\nCalculating densities...\n")
window_size <- 100000  # 100kb windows

calculate_density <- function(annotations, chr_info, window_size) {
    densities <- list()
    
    for (i in 1:nrow(chr_info)) {
        chr_name <- chr_info$chr[i]
        chr_len <- chr_info$length[i]
        
        # Create windows
        windows <- seq(1, chr_len, by = window_size)
        n_windows <- length(windows)
        
        # Get features on this chromosome
        chr_features <- annotations[annotations$chr == chr_name, ]
        
        # Count features in each window
        counts <- numeric(n_windows)
        for (j in 1:n_windows) {
            win_start <- windows[j]
            win_end <- min(win_start + window_size - 1, chr_len)
            
            # Count features overlapping this window
            overlapping <- chr_features$start <= win_end & chr_features$end >= win_start
            counts[j] <- sum(overlapping)
        }
        
        # Store results
        densities[[chr_name]] <- data.frame(
            chr = chr_name,
            start = windows,
            end = pmin(windows + window_size - 1, chr_len),
            count = counts,
            density = counts / (window_size / 1000)  # Features per kb
        )
    }
    
    return(do.call(rbind, densities))
}

gene_density <- calculate_density(genes, main_chrs, window_size)
te_density <- calculate_density(tes, main_chrs, window_size)

cat("Gene density range:", min(gene_density$density), "-", max(gene_density$density), "genes/kb\n")
cat("TE density range:", min(te_density$density), "-", max(te_density$density), "TEs/kb\n")

# Create circos plot
cat("\nCreating circos plot...\n")

pdf(output_pdf, width = 10, height = 10)

# Initialize circos
circos.clear()
circos.par(
    start.degree = 90,
    gap.degree = 2,
    track.margin = c(0.01, 0.01),
    cell.padding = c(0.02, 0, 0.02, 0)
)

# Set up chromosome layout
chr_df <- data.frame(
    chr = main_chrs$chr,
    start = 0,
    end = main_chrs$length
)

circos.initialize(factors = chr_df$chr, 
                  xlim = as.matrix(chr_df[, c("start", "end")]))

# Track 1: Chromosome ideogram
circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
        chr = CELL_META$sector.index
        xlim = CELL_META$xlim
        circos.text(mean(xlim), 0.5, chr, 
                   facing = "bending.inside", 
                   niceFacing = TRUE,
                   cex = 1.2, font = 2)
    },
    bg.col = "grey90",
    bg.border = NA,
    track.height = 0.08
)

# Track 2: Gene density
cat("Adding gene density track...\n")
circos.track(
    ylim = c(0, max(gene_density$density) * 1.1),
    panel.fun = function(x, y) {
        chr = CELL_META$sector.index
        chr_data = gene_density[gene_density$chr == chr, ]
        
        circos.lines(
            x = (chr_data$start + chr_data$end) / 2,
            y = chr_data$density,
            col = "#2E86AB",
            lwd = 1.5
        )
        
        # Add grid lines
        circos.yaxis(side = "left", labels.cex = 0.5)
    },
    bg.border = NA,
    track.height = 0.15
)

# Track 3: TE density
cat("Adding TE density track...\n")
circos.track(
    ylim = c(0, max(te_density$density) * 1.1),
    panel.fun = function(x, y) {
        chr = CELL_META$sector.index
        chr_data = te_density[te_density$chr == chr, ]
        
        circos.lines(
            x = (chr_data$start + chr_data$end) / 2,
            y = chr_data$density,
            col = "#F18F01",
            lwd = 1.5
        )
        
        # Add grid lines
        circos.yaxis(side = "left", labels.cex = 0.5)
    },
    bg.border = NA,
    track.height = 0.15
)

# Track 4: Gene annotations (as bars)
cat("Adding gene annotation track...\n")
circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
        chr = CELL_META$sector.index
        chr_genes = genes[genes$chr == chr, ]
        
        if (nrow(chr_genes) > 0) {
            # Subsample if too many genes
            if (nrow(chr_genes) > 1000) {
                chr_genes <- chr_genes[sample(nrow(chr_genes), 1000), ]
            }
            
            circos.rect(
                xleft = chr_genes$start,
                xright = chr_genes$end,
                ybottom = 0,
                ytop = 1,
                col = "#2E86AB",
                border = NA
            )
        }
    },
    bg.border = NA,
    track.height = 0.08
)

# Track 5: TE annotations (as bars)
cat("Adding TE annotation track...\n")
circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
        chr = CELL_META$sector.index
        chr_tes = tes[tes$chr == chr, ]
        
        if (nrow(chr_tes) > 0) {
            # Subsample if too many TEs
            if (nrow(chr_tes) > 1000) {
                chr_tes <- chr_tes[sample(nrow(chr_tes), 1000), ]
            }
            
            circos.rect(
                xleft = chr_tes$start,
                xright = chr_tes$end,
                ybottom = 0,
                ytop = 1,
                col = "#F18F01",
                border = NA
            )
        }
    },
    bg.border = NA,
    track.height = 0.08
)

# Add legend
legend(
    "center",
    legend = c("Gene Density", "TE Density", "Gene Annotations", "TE Annotations"),
    col = c("#2E86AB", "#F18F01", "#2E86AB", "#F18F01"),
    lty = c(1, 1, NA, NA),
    pch = c(NA, NA, 15, 15),
    lwd = c(2, 2, NA, NA),
    pt.cex = c(NA, NA, 2, 2),
    bty = "n",
    cex = 1.2
)

# Add title
title(main = "Gene and TE Annotation Landscape", cex.main = 1.5, font.main = 2)

circos.clear()
dev.off()

cat("\nPlot saved to:", output_pdf, "\n")

# Generate summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total genes:", nrow(genes), "\n")
cat("Total TEs:", nrow(tes), "\n")
cat("Genome size:", sum(main_chrs$length), "bp\n")
cat("Gene coverage:", sum(genes$end - genes$start + 1), "bp (", 
    round(sum(genes$end - genes$start + 1) / sum(main_chrs$length) * 100, 2), "%)\n")
cat("TE coverage:", sum(tes$end - tes$start + 1), "bp (",
    round(sum(tes$end - tes$start + 1) / sum(main_chrs$length) * 100, 2), "%)\n")

# Correlation between gene and TE density
cat("\n=== Density Correlation ===\n")
merged_density <- merge(gene_density, te_density, by = c("chr", "start", "end"))
correlation <- cor(merged_density$density.x, merged_density$density.y)
cat("Correlation between gene density and TE density:", round(correlation, 3), "\n")

if (correlation < -0.3) {
    cat("→ Negative correlation: TEs tend to be in gene-poor regions\n")
} else if (correlation > 0.3) {
    cat("→ Positive correlation: TEs and genes tend to co-occur\n")
} else {
    cat("→ Weak correlation: TEs and genes are independently distributed\n")
}

cat("\nDone!\n")