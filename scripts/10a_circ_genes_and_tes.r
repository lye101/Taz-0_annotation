# Load the circlize package
library(circlize)
library(tidyverse)
library(ComplexHeatmap)
suppressPackageStartupMessages(library(ComplexHeatmap)) 

setwd("/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/results")
# Load the TE annotation GFF3 file
gff_file <- "01a_EDTA_annotation/assembly.fasta.mod.EDTA.TEanno.gff3"
gff_data <- read.table(gff_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
genes_file <- "final/filtered.genes.renamed.gff3"
genes_gff <- read.table(genes_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
# Add "S" prefix to chromosome names to match FAI file
gff_data$V1 <- paste0("S", gff_data$V1)
genes_gff$V1 <- paste0("S", genes_gff$V1)

# Check the superfamilies present in the GFF3 file, and their counts
gff_data$V3 %>% table()
genes_gff$V3 %>% table()


# custom ideogram data
custom_ideogram <- read.table("/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/results/03a_fai_files/s.assembly.fasta.fai", header = FALSE, stringsAsFactors = FALSE)
custom_ideogram$chr <- custom_ideogram$V1
custom_ideogram$start <- 1
custom_ideogram$end <- custom_ideogram$V2
custom_ideogram <- custom_ideogram[, c("chr", "start", "end")]
custom_ideogram <- custom_ideogram[order(custom_ideogram$end, decreasing = T), ]
sum(custom_ideogram$end[1:20])

# Select only the first 20 longest scaffolds
custom_ideogram <- custom_ideogram[1:20, ]

# Keep only "gene" features
genes <- genes_gff %>%
    filter(V3 == "gene") %>%
    mutate(chrom = V1,
           start = V4,
           end = V5) %>%
    select(chrom, start, end) %>%
    filter(chrom %in% custom_ideogram$chr)

# Function to filter GFF3 data based on Superfamily
filter_superfamily <- function(gff_data, superfamily, custom_ideogram) {
    filtered_data <- gff_data %>%
        filter(V3 == superfamily, V1 %in% custom_ideogram$chr) %>%
        select(V1, V4, V5) %>%
        as.data.frame()
    
    # Return NULL if no data
    if (nrow(filtered_data) == 0) {
        message(paste("No data found for", superfamily))
        return(NULL)
    }
    
    # Ensure proper data types
    filtered_data[[1]] <- as.character(filtered_data[[1]])
    filtered_data[[2]] <- as.numeric(filtered_data[[2]])
    filtered_data[[3]] <- as.numeric(filtered_data[[3]])
    
    # Remove any rows with NA values or invalid coordinates
    filtered_data <- filtered_data[complete.cases(filtered_data), ]
    filtered_data <- filtered_data[filtered_data[[2]] < filtered_data[[3]], ]  # Ensure start < end
    
    if (nrow(filtered_data) == 0) {
        message(paste("No valid data after cleaning for", superfamily))
        return(NULL)
    }
    
    message(paste("Found", nrow(filtered_data), "features for", superfamily))
    
    return(filtered_data)
}

# Check data before plotting
gypsy_data <- filter_superfamily(gff_data, "Gypsy_LTR_retrotransposon", custom_ideogram)
copia_data <- filter_superfamily(gff_data, "Copia_LTR_retrotransposon", custom_ideogram)

hat_data <- filter_superfamily(gff_data, "hAT_TIR_transposon", custom_ideogram)
cacta_data <- filter_superfamily(gff_data, "CACTA_TIR_transposon", custom_ideogram)
mutator_data <- filter_superfamily(gff_data, "Mutator_TIR_transposon", custom_ideogram)
helitron_data <- filter_superfamily(gff_data, "helitron", custom_ideogram)

#Plot constants
track_height <- 0.1
window_size <- 1e6
countby <- "number"


pdf("TE_and_gene_density.pdf", width = 10, height = 10)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)

# Initialize the circos plot with the custom ideogram
circos.genomicInitialize(custom_ideogram)

# plot genes
circos.genomicDensity(genes,
                      count_by = countby,
                      col = "#3A7B00",
                      track.height = 0.13,
                      window.size = 1e5)

# Plot te density - only if data exists
if (!is.null(gypsy_data) && nrow(gypsy_data) > 0) {
    message("Plotting Gypsy data...")
    circos.genomicDensity(gypsy_data, col = c("#00800080"), track.height = track_height, window.size = window_size, count_by = countby)
}

if (!is.null(copia_data) && nrow(copia_data) > 0) {
    message("Plotting Copia data...")
    circos.genomicDensity(copia_data, col = c("#8B000080"), track.height = track_height, window.size = window_size, count_by = countby)
}

if (!is.null(hat_data) && nrow(hat_data) > 0) {
    circos.genomicDensity(hat_data, col = c("#FF000080"), track.height = track_height, window.size = window_size, count_by = countby)
}

if (!is.null(cacta_data) && nrow(cacta_data) > 0) {
    circos.genomicDensity(cacta_data, col = c("#0000FF80"), track.height = track_height, window.size = window_size, count_by = countby)
}

if (!is.null(mutator_data) && nrow(mutator_data) > 0) {
    circos.genomicDensity(mutator_data, col = c("#FFA50080"), track.height = track_height, window.size = window_size, count_by = countby)
}

if (!is.null(helitron_data) && nrow(helitron_data) > 0) {
    circos.genomicDensity(helitron_data, col = c("#80008080"), track.height = track_height, window.size = window_size, count_by = countby)
}

circos.clear()

lgd <- Legend(title = "TE Superfamilies \nand Genes", at = c("Genes","Gypsy", "Copia", "hAT", "CACTA", "Mutator", "Helitron"),
    legend_gp = gpar(fill = c("#3A7B00","#00800080", "#8B000080", "#FF000080", "#0000FF80", "#FFA50080", "#80008080"))
)

draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center", "center"))

dev.off()